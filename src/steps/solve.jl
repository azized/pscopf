using SplitApplyCombine

function tso_solve!(model_container::AbstractModelContainer,
                solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                dynamic_solving::Bool)
    if dynamic_solving
        @info "dynamic solve"
        @timeit TIMER_TRACKS "tso_dynamic_solve" iterative_solve_on_rso_constraints!(model_container,
                                                                            solve_fct, configs, uncertainties_at_ech, network)
    else
        @timeit TIMER_TRACKS "tso_modeling" begin
            @timeit TIMER_TRACKS "rso_cstrs" begin
                rso_combinations = available_combinations(network, get_target_timepoints(uncertainties_at_ech), get_scenarios(uncertainties_at_ech))
                @info @sprintf("adding %d constraints", length(rso_combinations))
                add_rso_flows_exprs!(model_container,
                                    rso_combinations,
                                    uncertainties_at_ech,
                                    network)
                add_rso_constraints!(model_container,
                                    rso_combinations,
                                    network)
            end
        end
        @info "tso_solve!: actual solve"
        @timeit TIMER_TRACKS "tso_solve" solve_fct(model_container, configs)
    end
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
    consumptions = build_fast_access_consumptions(uncertainties_at_ech, network)
    ptdf_l = build_fast_ptdf(network.ptdf)

    did_add_constraints = true
    while(did_add_constraints)
        @info @sprintf("dynamic solve: handling %d constraints out of %d\n",
            length(get_rso_constraints(model_container)), theoretical_nb_combinations(network, timepoints, scenarios))
        solve_fct(model_container, configs)

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
                                                                                                consumptions, ptdf_l)
    end

    @info @sprintf("Dynamically added %d constraints out of %d possible combinations\n",
            length(get_rso_constraints(model_container)), theoretical_nb_combinations(network, timepoints, scenarios))
end


function generate_rso_constraints!(model_container::AbstractModelContainer,
                                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                                timepoints, scenarios,
                                max_add_per_iter,
                                consumptions, fast_ptdf)

    violated_combinations = compute_violated_combinations(model_container, consumptions, network, timepoints, scenarios, fast_ptdf)
    has_violations = !isempty(violated_combinations)

    @timeit TIMER_TRACKS "sort_violations" violated_combinations_to_add = violations_to_add(violated_combinations, max_add_per_iter)

    # set start values for the next iteration before invalidating the model
    if !isempty(violated_combinations_to_add)
        set_start_values!(get_model(model_container))
    end

    timed_l = @timed add_constraints!(model_container,
                        violated_combinations_to_add, uncertainties_at_ech, network)
    @info @sprintf("adding RSO constraints took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))

    return has_violations
end


function compute_violated_combinations(model_container::AbstractModelContainer,
                                    uncertainties_at_ech, network,
                                    timepoints, scenarios,
                                    fast_ptdf
                                    )::Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}
    violated_combinations = Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}()
    timed_l = @timed  begin
        p_values_lim = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_p_injected(model_container, Networks.LIMITABLE))], keys(get_p_injected(model_container, Networks.LIMITABLE))))
        p_values_pil = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_p_injected(model_container, Networks.PILOTABLE))], keys(get_p_injected(model_container, Networks.PILOTABLE))))
        p_values_lol = value.(JuMP.Containers.DenseAxisArray([v for v in values(get_local_lol(model_container))], keys(get_local_lol(model_container))))
    end
    @info @sprintf("creating values container took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))
    timed_l = @timed for (branch_id,ts,s,network_case) in available_combinations(network, timepoints, scenarios)
        # constraint not considered in the model yet
        if !((branch_id,ts,s,network_case) in keys(get_rso_constraints(model_container)))
            branch = Networks.get_branch(network, branch_id)
            @timeit TIMER_TRACKS "create_or_get_expr" flow_val_l = flow_val(branch, ts, s, network_case,
                                                                        uncertainties_at_ech, network,
                                                                        p_values_lim, p_values_pil, p_values_lol, fast_ptdf)

            branch_limit = Networks.safeget_limit(branch, network_case)
            @timeit TIMER_TRACKS "eval_expr" violation_l = max(0., abs(flow_val_l) - branch_limit)
            if violation_l > 0
                # store violated combinations to add constraints later, not to invalidate the model now
                push!(violated_combinations, (branch, ts, s, network_case) => violation_l)
            end
        end
    end
    @info @sprintf("verifying RSO constraint violations took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))
    @info @sprintf("number of violated constraints : %d", length(violated_combinations))

    return violated_combinations
end

function violations_to_add_by_ts(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    #sort by violation #Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}
    violated_combinations_to_add = Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}()
    for grp in group( pair_combination_violation -> pair_combination_violation[1][2], violated_combinations)
        sorted_grp_violations = sort(grp, by=pair_combination_violation->pair_combination_violation[2], rev=true)
        append!(violated_combinations_to_add, sorted_grp_violations[1:min(max_add_per_iter,end)])
    end
    return violated_combinations_to_add
end

function violations_to_add(violated_combinations, max_add_per_iter)::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}}
    #sort by violation
    sorted_violations = sort(collect(violated_combinations), by= x -> x[2], rev=true)

    violated_combinations_to_add = sorted_violations[1: min(max_add_per_iter, length(sorted_violations))]
    return violated_combinations_to_add
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
