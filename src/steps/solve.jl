function tso_solve!(model_container::AbstractModelContainer,
                solve_fct::Base.Callable, configs::AbstractRunnableConfigs,
                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network,
                dynamic_solving::Bool)
    if dynamic_solving
        @timeit TIMER_TRACKS "tso_dynamic_solve" iterative_solve_on_rso_constraints!(model_container,
                                                                            solve_fct, configs, uncertainties_at_ech, network)
    else
        @debug "adding RSO constraints"
        @timeit TIMER_TRACKS "tso_modeling" begin
            add_rso_flows_exprs!(model_container,
                                get_rso_combinations(model_container),
                                uncertainties_at_ech,
                                network)
            add_rso_constraints!(model_container,
                                get_rso_combinations(model_container),
                                network)
        end
        @debug "tso_solve!: actual solve"
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

        did_add_constraints = generate_rso_constraints!(model_container, uncertainties_at_ech, network)
    end

    @info @sprintf("Dynamically added %d constraints out of %d possible combinations\n",
            length(get_rso_constraints(model_container)), length(get_rso_combinations(model_container)))
end


function generate_rso_constraints!(model_container::AbstractModelContainer,
                                uncertainties_at_ech::UncertaintiesAtEch, network::Networks.Network)
    max_add_per_iter = get_config("MAX_ADD_RSO_CSTR_PER_ITER")
    violated_combinations_to_add = Vector{Tuple{Networks.Branch,DateTime,String,String}}()
    for (branch_id,ts,s,network_case) in get_rso_combinations(model_container)
        # constraint not considered in the model yet
        if !((branch_id,ts,s,network_case) in keys(get_rso_constraints(model_container)))
            branch = Networks.get_branch(network, branch_id)
            flow_expr_l = flow_expr(model_container,
                                    branch, ts, s, network_case,
                                    uncertainties_at_ech, network)

            branch_limit = Networks.safeget_limit(branch, network_case)
            if abs(value(flow_expr_l)) > branch_limit
                add_rso_flow_expr!(model_container,
                                flow_expr_l,
                                branch, ts, s, network_case)
                # store violated combinations to add constraints later, not to invalidate the model now
                push!(violated_combinations_to_add, (branch, ts, s, network_case))
            end
        end

        if !isnothing(max_add_per_iter) && (length(violated_combinations_to_add) >= max_add_per_iter)
            break;
        end
    end

    if !isempty(violated_combinations_to_add)
        set_start_values!(get_model(model_container))
    end

    for (branch, ts, s, network_case) in violated_combinations_to_add
        add_rso_constraint!(get_model(model_container), get_rso_constraints(model_container),
                        get_flows(model_container),
                        branch, ts, s, network_case)
    end

    has_violations = !isempty(violated_combinations_to_add)
    return has_violations
end
