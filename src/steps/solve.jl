using SplitApplyCombine

function tso_solve!(model_container::AbstractModelContainer,
                solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                dynamic_solving::Bool)
    full_solve_timed_l = @timed begin
        if dynamic_solving
            @info "dynamic solve"
            @timeit TIMER_TRACKS "tso_dynamic_solve" iterative_solve_on_rso_constraints!(model_container,
                                                                                solve_fct, configs, uncertainties_at_ech, network)
        else
            @timeit TIMER_TRACKS "tso_modeling" begin
                @timeit TIMER_TRACKS "rso_cstrs" begin
                    @info @sprintf("adding %d constraints", sum(1 for iter in available_combinations(network, get_target_timepoints(uncertainties_at_ech), get_scenarios(uncertainties_at_ech))))
                    add_rso_flows_exprs!(model_container,
                                        available_combinations(network, get_target_timepoints(uncertainties_at_ech), get_scenarios(uncertainties_at_ech)),
                                        uncertainties_at_ech,
                                        network)
                    add_rso_constraints!(model_container,
                                        available_combinations(network, get_target_timepoints(uncertainties_at_ech), get_scenarios(uncertainties_at_ech)),
                                        network)
                end
            end
            @info "tso_solve!: actual solve"
            @timeit TIMER_TRACKS "tso_solve" solve_fct(model_container, configs)
            TSO_SOLVE_RECORDS.nb_iters = 1
        end
    end

    #TODO : avoid hardcode
    TSO_SOLVE_RECORDS.usecase = substring_from(filter(!isspace, configs.out_path), "pscopf")
    TSO_SOLVE_RECORDS.dynamic = dynamic_solving
    TSO_SOLVE_RECORDS.nb_added_constraints = length(get_rso_constraints(model_container))
    TSO_SOLVE_RECORDS.nb_total_constraints = theoretical_nb_combinations(network, get_target_timepoints(uncertainties_at_ech), get_scenarios(uncertainties_at_ech))
    TSO_SOLVE_RECORDS.total_cstr_generation_time = TimerOutputs.ncalls(dynamic_solving ?
                                                                    TIMER_TRACKS["run_model"]["tso_dynamic_solve"]["generate_rso_cstrs"] :
                                                                    TIMER_TRACKS["run_model"]["tso_modeling"]["rso_cstrs"] )
    TSO_SOLVE_RECORDS.total_solving_time = full_solve_timed_l.time
    TSO_SOLVE_RECORDS.total_cstr_generation_time = TimerOutputs.time(dynamic_solving ?
                                                                    TIMER_TRACKS["run_model"]["tso_dynamic_solve"]["generate_rso_cstrs"] :
                                                                    TIMER_TRACKS["run_model"]["tso_modeling"]["rso_cstrs"]) /1e9
    write_record!(TSO_SOLVE_RECORDS)
end

function build_fast_ptdf(ptdf::PTDFDict)
    result_ptdf = Dict{Tuple{String, String, String},Float64}()
    for (case, ptdf_vals) in ptdf
        for (branch, ptdf_elts) in ptdf_vals
            for (bus, val) in ptdf_elts
                result_ptdf[case, branch, bus] = val
            end
        end
    end
    return result_ptdf
end

function build_fast_branches(branches::SortedDict{String, Networks.Branch})
    result_limits = Dict{Tuple{String, String}, Tuple{Networks.Branch, Float64}}()
    for (branch_id, branch) in branches
        for (network_case, val_limit) in get_limit(branch)
            result_limits[branch_id, network_case] = (branch, val_limit)
        end
    end
    return result_limits
end

function build_fast_access_consumptions(uncertainties_at_ech, network)
    result = Dict{Tuple{String, DateTime, String}, Float64}()
    for (injection_id, injection_uncertainties) in uncertainties_at_ech
        if !haskey(network.buses, injection_id)
            continue # skip limitable injections
        end

        for (ts, by_s_uncertainties) in injection_uncertainties
            for (s, val) in by_s_uncertainties
                result[injection_id, ts, s] = val
            end
        end
    end
    return result
end

function get_uncertainties(uncertainties, injection_id, ts, s)
    return uncertainties[injection_id, ts, s]
end


function iterative_solve_on_rso_constraints!(model_container::AbstractModelContainer,
                                        solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                                        uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network)
    timepoints = get_target_timepoints(uncertainties_at_ech)
    scenarios = get_scenarios(uncertainties_at_ech)
    consumptions = build_fast_access_consumptions(uncertainties_at_ech, network) #TODO : compute once directly as base flows
    ptdf_l = build_fast_ptdf(network.ptdf) #TODO : once early on
    fast_branches_l = build_fast_branches(network.branches) #TODO : once early on

    did_add_constraints = true
    while(did_add_constraints)
        @info @sprintf("dynamic solve: handling %d constraints out of %d\n",
            length(get_rso_constraints(model_container)), theoretical_nb_combinations(network, timepoints, scenarios))
        solver_timed_l = @timed solve_fct(model_container, configs)

        # model is infeasible, adding constraints will not solve the issue
        if get_status(model_container) in [pscopf_INFEASIBLE, pscopf_HAS_SLACK]
            @warn("check status : may be due to numerical errors")
            return
        end

        max_add_per_iter = get_config("MAX_ADD_RSO_CSTR_PER_ITER") < 0 ? theoretical_nb_combinations(network, timepoints, scenarios) : get_config("MAX_ADD_RSO_CSTR_PER_ITER")
        @timeit TIMER_TRACKS "generate_rso_cstrs" did_add_constraints = generate_rso_constraints!(model_container,
                                                                                                uncertainties_at_ech, network,
                                                                                                timepoints, scenarios,
                                                                                                max_add_per_iter,
                                                                                                consumptions, ptdf_l, fast_branches_l)

        DYNAMIC_SOLVE_RECORDS.nb_added_constraints = length(get_rso_constraints(model_container))
        DYNAMIC_SOLVE_RECORDS.nb_total_constraints = theoretical_nb_combinations(network, timepoints, scenarios)
        DYNAMIC_SOLVE_RECORDS.solver_time = solver_timed_l.time
        write_record!(DYNAMIC_SOLVE_RECORDS)
    end

    TSO_SOLVE_RECORDS.nb_iters = DYNAMIC_SOLVE_RECORDS.iter
    @info @sprintf("Dynamically added %d constraints out of %d possible combinations\n",
            length(get_rso_constraints(model_container)), theoretical_nb_combinations(network, timepoints, scenarios))
end


function generate_rso_constraints!(model_container::AbstractModelContainer,
                                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                                timepoints, scenarios,
                                max_add_per_iter,
                                consumptions, fast_ptdf, fast_branches_l)

    @timeit TIMER_TRACKS "compute_violations" begin
        timed_l = @timed violated_combinations = compute_violated_combinations(model_container, consumptions, network, timepoints, scenarios, fast_ptdf, fast_branches_l)
        @info @sprintf("compute_violations took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))
        has_violations = !isempty(violated_combinations)

        DYNAMIC_SOLVE_RECORDS.nb_violated_constraints = length(violated_combinations)
    end

    # set start values for the next iteration before invalidating the model
    if !isempty(violated_combinations)
        set_start_values!(get_model(model_container))
    end

    @timeit TIMER_TRACKS "sort_violations" begin
        violations_filtering_time_l = @timed violated_combinations_to_add = violations_to_add(violated_combinations, max_add_per_iter)
    end

    timed_l = @timed add_constraints!(model_container,
                        violated_combinations_to_add, uncertainties_at_ech, network)
    @info @sprintf("adding RSO constraints took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))

    DYNAMIC_SOLVE_RECORDS.violations_filtering_time = violations_filtering_time_l.time

    return has_violations
end

function compute_p_values_by_bus(model_container,
                                uncertainties_at_ech, network, timepoints, scenarios
                                )::Dict{String, Dict{Tuple{DateTime,String}, Float64}}
    p_values = Dict{String, Dict{Tuple{DateTime,String}, Float64}}()

    timed_l = @timed  begin
        p_values_lim = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_p_injected(model_container, Networks.LIMITABLE))], keys(get_p_injected(model_container, Networks.LIMITABLE))))
        p_values_pil = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_p_injected(model_container, Networks.PILOTABLE))], keys(get_p_injected(model_container, Networks.PILOTABLE))))
        p_values_lol = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_local_lol(model_container))], keys(get_local_lol(model_container))))
    end
    @info @sprintf("creating values container took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))

    timed_l = @timed for (bus_id, bus) in network.buses
        bus_p_values = get!(p_values, bus_id, Dict{Tuple{DateTime,String}, Float64}())

        for ts in timepoints, s in scenarios
            bus_p_values[ts,s] = - get_uncertainties(uncertainties_at_ech, bus_id, ts, s)
            bus_p_values[ts,s] += p_values_lol[(bus_id, ts, s)]
        end

        for gen in Networks.get_generators(bus)
            gen_id = get_id(gen)
            gen_type = get_type(gen)

            if gen_type == LIMITABLE
                for ts in timepoints, s in scenarios
                    bus_p_values[ts,s] += p_values_lim[(gen_id, ts, s)]
                end
            end

            if gen_type == PILOTABLE
                for ts in timepoints, s in scenarios
                    bus_p_values[ts,s] += p_values_pil[(gen_id, ts, s)]
                end
            end
        end

    end
    @info @sprintf("creating p_values by bus took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))

    return p_values
end

function compute_flows(model_container,
                    uncertainties_at_ech, network,
                    timepoints, scenarios,
                    fast_ptdf)

    p_values = compute_p_values_by_bus(model_container, uncertainties_at_ech, network, timepoints, scenarios)

    flow_computation_time_l = @timed begin
        flows_l = Dict{Tuple{String,DateTime,String,String}, Float64}(combination => 0. for combination in available_combinations(network, timepoints, scenarios))
        for ((network_case, branch_id, bus_id), ptdf_val_l) in fast_ptdf
            for ((ts,s), p_val) in p_values[bus_id]
                # if !((branch_id,ts,s,network_case) in keys(get_rso_constraints(model_container)))
                flows_l[branch_id,ts,s,network_case] += p_val * ptdf_val_l
            end
        end
    end
    @info @sprintf("computing flows took %s s and allocated %s kB", flow_computation_time_l.time, (flow_computation_time_l.bytes/1024))

    DYNAMIC_SOLVE_RECORDS.flow_computation_time = flow_computation_time_l.time

    return flows_l
end

function compute_violated_combinations(model_container::AbstractModelContainer,
                                    uncertainties_at_ech, network,
                                    timepoints, scenarios,
                                    fast_ptdf, fast_branches_l
                                    )::Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}
    violated_combinations = Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}()

    flows_l = compute_flows(model_container, uncertainties_at_ech, network, timepoints, scenarios, fast_ptdf)

    flow_verification_time_l = @timed for ((branch_id,ts,s,network_case),flow_val) in flows_l
        branch_limit = fast_branches_l[branch_id, network_case][2]
        violation_l = max(0., abs(flow_val) - branch_limit)
        if violation_l > 1e-09
            # store violated combinations to add constraints later, not to invalidate the model now
            push!(violated_combinations, (fast_branches_l[branch_id, network_case][1], ts, s, network_case) => violation_l)
        end
    end
    @info @sprintf("verifying RSO constraint violations took %s s and allocated %s kB", flow_verification_time_l.time, (flow_verification_time_l.bytes/1024))

    @info @sprintf("number of violated constraints : %d", length(violated_combinations))

    DYNAMIC_SOLVE_RECORDS.flow_verification_time = flow_verification_time_l.time

    return violated_combinations
end

function violations_to_add_by_ts_group(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    println(max_add_per_iter)
    #max_violations[ts][branch,s] => (violation, network_case)
    dict_max_violations = Dict{DateTime, Dict{Tuple{Networks.Branch,String},Tuple{Float64,String}}}()

    for ((br,ts,s,nc), val) in violated_combinations
        dict_max_violations_at_ts = get!(dict_max_violations, ts,
                                        Dict{Tuple{Networks.Branch,String},Tuple{Float64,String}}())
        if val > get(dict_max_violations_at_ts, (br,s), (0.,nothing))[1]
            dict_max_violations_at_ts[br,s] = (val, nc)
        end
    end

    dict_sorted_max_violations = Dict{DateTime, Vector{Pair{Tuple{Networks.Branch,String},Tuple{Float64,String}}}}()
    for (ts, dict_max_violations_at_ts) in dict_max_violations
        println(ts)
        dict_sorted_max_violations[ts] = sort(collect(dict_max_violations_at_ts),
                                            by=x->x[2][1], rev=true)
    end

    violated_combinations_to_add = Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}()
    for (ts, sorted_violations_at_ts) in dict_sorted_max_violations
        for ((branch,s),(val, case)) in sorted_violations_at_ts[1:min(end, max_add_per_iter)]
            push!(violated_combinations_to_add, (branch,ts,s,case) => val)
        end
    end

    return violated_combinations_to_add
end

function violations_to_add_by_ts(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    #sort by violation #Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}
    violated_combinations_to_add = Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}()
    by_ts = pair_combination_violation -> pair_combination_violation[1][2]
    by_violation = pair_combination_violation -> pair_combination_violation[2]
    for grp in group( by_ts, violated_combinations)
        sorted_grp_violations = sort(grp, by=by_violation, rev=true)
        append!(violated_combinations_to_add, sorted_grp_violations[1:min(max_add_per_iter,end)])
    end
    return violated_combinations_to_add
end

function violations_to_add_by_highest(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    #sort by violation
    sorted_violations = sort(collect(violated_combinations), by= x -> x[2], rev=true)

    violated_combinations_to_add = sorted_violations[1: min(max_add_per_iter, length(sorted_violations))]
    return violated_combinations_to_add
end

function violations_to_add(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    if (get_config("VIOLATIONS_ADD_FUNCTION") == "HIGHEST")
        return violations_to_add_by_highest(violated_combinations, max_add_per_iter)
    elseif (get_config("VIOLATIONS_ADD_FUNCTION") == "TS")
        return violations_to_add_by_ts(violated_combinations, max_add_per_iter)
    elseif (get_config("VIOLATIONS_ADD_FUNCTION") == "TS_GROUP")
        return violations_to_add_by_ts_group(violated_combinations, max_add_per_iter)
    else
        error("unrecognized violations adding function "*get_config("VIOLATIONS_ADD_FUNCTION"))
    end
end

function add_constraints!(model_container::AbstractModelContainer,
                        violated_combinations_to_add,
                        uncertainties_at_ech, network)
    @timeit TIMER_TRACKS "add_cstrs" for ((branch, ts, s, network_case) ,_) in violated_combinations_to_add
        flow_expr_l = flow_expr(model_container,
                                branch, ts, s, network_case,
                                uncertainties_at_ech, network)
        add_rso_flow_expr!(model_container,
                                flow_expr_l,
                                branch, ts, s, network_case)
        add_rso_constraint!(get_model(model_container), get_rso_constraints(model_container),
                        get_flows(model_container),
                        branch, ts, s, network_case)
    end
end
