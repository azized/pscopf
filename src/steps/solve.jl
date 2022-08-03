function tso_solve!(model_container::AbstractModelContainer,
                solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                dynamic_solving::Bool)
    if dynamic_solving
        @info "dynamic solve"
        @timeit TIMER_TRACKS "tso_dynamic_solve" iterative_solve_on_rso_constraints!(model_container,
                                                                            solve_fct, configs, uncertainties_at_ech, network)
    else
        @info @sprintf("adding %d RSO constraints", length(get_rso_combinations(model_container)))
        @timeit TIMER_TRACKS "tso_modeling" begin
        @timeit TIMER_TRACKS "rso_cstrs" begin
            add_rso_flows_exprs!(model_container,
                                get_rso_combinations(model_container),
                                uncertainties_at_ech,
                                network)
            add_rso_constraints!(model_container,
                                get_rso_combinations(model_container),
                                network)
        end
        end
        @info "tso_solve!: actual solve"
        @timeit TIMER_TRACKS "tso_solve" solve_fct(model_container, configs)
    end
end

function iterative_solve_on_rso_constraints!(model_container::AbstractModelContainer,
                                        solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                                        uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network)
    did_add_constraints = true
    while(did_add_constraints)
        @info @sprintf("dynamic solve: handling %d constraints out of %d\n",
            length(get_rso_constraints(model_container)), length(get_rso_combinations(model_container)))
        solve_fct(model_container, configs)

        # model is infeasible, adding constraints will not solve the issue
        if get_status(model_container) in [pscopf_INFEASIBLE, pscopf_HAS_SLACK]
            @warn("check status : may be due to numerical errors")
            return
        end

        @timeit TIMER_TRACKS "generate_rso_cstrs" did_add_constraints = generate_rso_constraints!(model_container, uncertainties_at_ech, network)
    end

    @info @sprintf("Dynamically added %d constraints out of %d possible combinations\n",
            length(get_rso_constraints(model_container)), length(get_rso_combinations(model_container)))
end


function generate_rso_constraints!(model_container::AbstractModelContainer,
                                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network)
    max_add_per_iter = get_config("MAX_ADD_RSO_CSTR_PER_ITER") < 0 ? length(get_rso_combinations(model_container)) : get_config("MAX_ADD_RSO_CSTR_PER_ITER")

    violated_combinations = Dict{Tuple{Networks.Branch,DateTime,String,String}, Float64}()
    timed_l = @timed for (branch_id,ts,s,network_case) in get_rso_combinations(model_container)
        # constraint not considered in the model yet
        if !((branch_id,ts,s,network_case) in keys(get_rso_constraints(model_container)))
            branch = Networks.get_branch(network, branch_id)
            @timeit TIMER_TRACKS "create_or_get_expr" flow_expr_l = flow_expr(model_container,
                                    branch, ts, s, network_case,
                                    uncertainties_at_ech, network)

            branch_limit = Networks.safeget_limit(branch, network_case)
            @timeit TIMER_TRACKS "eval_expr" violation_l = max(0., abs(value(flow_expr_l)) - branch_limit)
            if violation_l > 0
                # store violated combinations to add constraints later, not to invalidate the model now
                push!(violated_combinations, (branch, ts, s, network_case) => violation_l)
            end
        end
    end
    @info @sprintf("verifying RSO constraint violations took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))
    @info @sprintf("number of violated constraints : %d", length(violated_combinations))
    has_violations = !isempty(violated_combinations)

    @timeit TIMER_TRACKS "sort_violations" ( sorted_violations::Vector{Pair{Tuple{Networks.Branch,DateTime,String,String}, Float64}} =
                        sort(collect(violated_combinations), by= x -> x[2], rev=true) )
    violated_combinations_to_add = sorted_violations[1: min(max_add_per_iter, length(sorted_violations))]


    # set start values for the next iteration before invalidating the model
    if !isempty(violated_combinations_to_add)
        set_start_values!(get_model(model_container))
    end

    timed_l = @timed begin
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
    @info @sprintf("adding RSO constraints took %s s and allocated %s kB", timed_l.time, (timed_l.bytes/1024))

    return has_violations
end
