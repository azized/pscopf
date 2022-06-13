using .Networks

using Dates
using DataStructures
using JuMP
using Printf
using Parameters

@with_kw mutable struct RSOAssessmentModel <: AbstractModelContainer
    model::Model = Model()

    # Uncertainty model
    #bus,ts
    uncertain_load::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()
    #limitableUnit,ts
    uncertain_prod::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()

    # Market model
    #unit, ts
    p_injected::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()
    #pilotable, ts
    b_in::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()
    b_marg::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()
    b_out::SortedDict{Tuple{String, DateTime},VariableRef} =
        SortedDict{Tuple{String, DateTime},VariableRef}()

    #ts : rho (subproblem variable)
    overflow::SortedDict{DateTime,VariableRef} =
        SortedDict{DateTime,VariableRef}()

    # duals : branch,ts (subproblem)
    pos_flow_limit_duals::SortedDict{Tuple{String,DateTime},VariableRef} =
        SortedDict{Tuple{String,DateTime},VariableRef}()
    pos_flow_limit_indicators::SortedDict{Tuple{String,DateTime},VariableRef} =
        SortedDict{Tuple{String,DateTime},VariableRef}()
    neg_flow_limit_duals::SortedDict{Tuple{String,DateTime},VariableRef} =
        SortedDict{Tuple{String,DateTime},VariableRef}()
    neg_flow_limit_indicators::SortedDict{Tuple{String,DateTime},VariableRef} =
        SortedDict{Tuple{String,DateTime},VariableRef}()

    #sum_bus ptdf[branch,bus] * (prod_bus - conso_bus)
    signed_flow::SortedDict{Tuple{String,DateTime},AffExpr} =
        SortedDict{Tuple{String,DateTime},AffExpr}()

    #objective
    full_obj::AffExpr = AffExpr(0.)
end

function add_load_uncertainties_vars!(model_container::RSOAssessmentModel,
                                    network, TS, assessment_uncertainties)
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        lb_l = get_assessment_uncertainties_lb(assessment_uncertainties, bus_id)
        ub_l = get_assessment_uncertainties_ub(assessment_uncertainties, bus_id)
        for ts in TS
            name =  @sprintf("uncertain_load[%s,%s]", bus_id, ts);
            model_container.uncertain_load[bus_id,ts] = @variable(model_container.model, base_name = name,
                                                                lower_bound = lb_l, upper_bound = ub_l)
        end
    end

    return model_container.uncertain_load
end
function add_prod_uncertainties_vars!(model_container::RSOAssessmentModel,
                                    network, TS, assessment_uncertainties, tso_actions)
    for limitable_gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(limitable_gen)
        lb_l = get_assessment_uncertainties_lb(assessment_uncertainties, gen_id)
        for ts in TS
            # this allows to set p_inj[gen_id] == uncertain_prod without needing to express min(uncertain, limit) in the model
            ub_l = min(get_assessment_uncertainties_ub(assessment_uncertainties, gen_id),
                    get_limitation(tso_actions, gen_id, ts),
                    Networks.get_p_max(limitable_gen))
            name =  @sprintf("uncertain_prod[%s,%s]", gen_id, ts);
            model_container.uncertain_prod[gen_id, ts] = @variable(model_container.model, base_name=name,
                                                                lower_bound = lb_l, upper_bound = ub_l)
        end
    end

    return model_container.uncertain_prod
end
function add_uncertainties_vars!(model_container::RSOAssessmentModel,
                                network, TS, assessment_uncertainties, tso_actions)
    add_load_uncertainties_vars!(model_container,
                                network, TS, assessment_uncertainties)
    add_prod_uncertainties_vars!(model_container,
                                network, TS, assessment_uncertainties, tso_actions)
    return model_container
end


function add_limitable_prod_vars!(model_container::RSOAssessmentModel, network, TS)
    for limitable_gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(limitable_gen)
        p_max = Networks.get_p_max(limitable_gen)
        for ts in TS
            name =  @sprintf("p_injected[%s,%s]", gen_id, ts)
            model_container.p_injected[gen_id,ts] =
                @variable(model_container.model, base_name=name, lower_bound=0., upper_bound=p_max)
        end
    end

    return model_container.p_injected
end

function add_use_limitables_constraints!(model_container::RSOAssessmentModel, network, TS, tso_actions)
    for limitable_gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(limitable_gen)
        for ts in TS
            limit_l = get_limitation(tso_actions, gen_id, ts)
            @assert(upper_bound(model_container.uncertain_prod[gen_id, ts]) <= limit_l)
            @constraint(model_container.model,
                        model_container.p_injected[gen_id,ts] == model_container.uncertain_prod[gen_id, ts])
        end
    end
    return model_container
end

function add_pilotable_prod_vars!(model_container_p::RSOAssessmentModel, network, TS, tso_actions)
    p_injected_l = model_container_p.p_injected
    b_in_l = model_container_p.b_in
    b_marg_l = model_container_p.b_marg
    b_out_l = model_container_p.b_out

    for pilotable_gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        gen_id = Networks.get_id(pilotable_gen)
        for ts in TS
            pmin_l, pmax_l = safeget_imposition(tso_actions, gen_id, ts)
            name_l =  @sprintf("b_in[%s,%s]", gen_id, ts);
            b_in_l[gen_id,ts] = @variable(model_container_p.model, base_name = name_l, binary = true)
            name_l =  @sprintf("b_marginal[%s,%s]", gen_id, ts);
            b_marg_l[gen_id,ts] = @variable(model_container_p.model, base_name = name_l, binary = true)
            name_l =  @sprintf("b_out[%s,%s]", gen_id, ts);
            b_out_l[gen_id,ts] = @variable(model_container_p.model, base_name = name_l, binary = true)
            name_l =  @sprintf("p_injected[%s,%s]", gen_id, ts);
            p_injected_l[gen_id,ts] = @variable(model_container_p.model, base_name = name_l,
                                                lower_bound=0., upper_bound = pmax_l);

            @constraint(model_container_p.model,
                        p_injected_l[gen_id,ts] >= ( pmax_l * b_in_l[gen_id,ts] +  pmin_l * b_marg_l[gen_id,ts] ));
            @constraint(model_container_p.model,
                        p_injected_l[gen_id,ts] <= (1 - b_out_l[gen_id,ts]) * pmax_l);
            @constraint(model_container_p.model,
                        b_out_l[gen_id,ts] + b_marg_l[gen_id,ts] + b_in_l[gen_id,ts] == 1.);
        end
    end

    return model_container_p
end

function add_market_vars!(model_container::RSOAssessmentModel, network, TS, tso_actions)
    add_limitable_prod_vars!(model_container, network, TS)
    add_pilotable_prod_vars!(model_container, network, TS, tso_actions)
end

function add_cheapest_prod_constraints!(model_container::RSOAssessmentModel, network, TS)
    pilotables = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    for (index, gen_1) in enumerate(pilotables)
        cost_1 = Networks.get_prop_cost(gen_1)
        gen_id_1 = Networks.get_id(gen_1)
        for gen_2 in pilotables[index:end]
            cost_2 = Networks.get_prop_cost(gen_2)
            if cost_1 < cost_2
                gen_id_2 = Networks.get_id(gen_2)
                for ts in TS
                    @constraint(model_container.model, model_container.b_in[gen_id_1,ts] >= model_container.b_in[gen_id_2,ts] );
                end
            end
        end
    end

    return model_container
end

function add_eod_constraint!(model_container_p::RSOAssessmentModel, network, TS)

    for ts in TS

        supply_l = AffExpr(0)
        for gen in Networks.get_generators(network)
            gen_id = Networks.get_id(gen)
            supply_l += model_container_p.p_injected[gen_id, ts];
        end

        demand_l = AffExpr(0.)
        for bus in Networks.get_buses(network)
            bus_id = Networks.get_id(bus)
            demand_l += model_container_p.uncertain_load[bus_id, ts];
        end

        @constraint(model_container_p.model, supply_l == demand_l);
    end

    return model_container_p
end

function add_overflow_var!(model_container::RSOAssessmentModel, TS, upper_bound)
    for ts in TS
        name = @sprintf("overflow[%s]", ts)
        model_container.overflow[ts] =
            @variable(model_container.model, base_name=name, lower_bound=0., upper_bound=upper_bound)
    end
    return model_container.overflow
end

function create_objective!(model_container::RSOAssessmentModel, network, TS)::AffExpr
    model_container.full_obj = AffExpr(0.)
    add_coeffxsum_cost!(model_container.full_obj, model_container.overflow, 1.)
    return model_container.full_obj
end

function add_uncertainties_constraints!(model_container_l::RSOAssessmentModel, network, TS, assessment_uncertainties, configs)
    #FIXME constraints on uncertainties go here !
    #e.g | u[bus1] - u[bus2] | < 5
    #e.g ( 1 - u[bus1] / ub[limitable_1] ) <  ( 1 - u[limitable_1] / ub[limitable_1] )
end

function create_flow_limit_duals_and_indicators(model_container::RSOAssessmentModel, network, TS)
    for branch in Networks.get_branches(network)
        branch_id = Networks.get_id(branch)
        for ts in TS
            name=@sprintf("pos_flow_duals[%s,%s]", branch_id, ts)
            model_container.pos_flow_limit_duals[branch_id, ts] =
                @variable(model_container.model, base_name=name, lower_bound=0., upper_bound=1.)
            name=@sprintf("pos_flow_indicators[%s,%s]", branch_id, ts)
            model_container.pos_flow_limit_indicators[branch_id, ts] =
                @variable(model_container.model, base_name=name, binary=true)
            name=@sprintf("neg_flow_duals[%s,%s]", branch_id, ts)
            model_container.neg_flow_limit_duals[branch_id, ts] =
                @variable(model_container.model, base_name=name, lower_bound=0., upper_bound=1.)
            name=@sprintf("neg_flow_indicators[%s,%s]", branch_id, ts)
            model_container.neg_flow_limit_indicators[branch_id, ts] =
                @variable(model_container.model, base_name=name, binary=true)
        end
    end

    return model_container
end

function create_flow_expressions(model_container::RSOAssessmentModel, network, TS)
    for branch in Networks.get_branches(network)
        branch_id = Networks.get_id(branch)

        for ts in TS

            branch_flow_expr = AffExpr(0.)
            for bus in Networks.get_buses(network)
                bus_id = Networks.get_id(bus)
                ptdf = Networks.safeget_ptdf_elt(network, branch_id, bus_id)

                # + injections
                for gen in Networks.get_generators(bus)
                    gen_id = Networks.get_id(gen)
                    branch_flow_expr += ptdf * model_container.p_injected[gen_id,ts]
                end

                # - loads
                branch_flow_expr -= ptdf * model_container.uncertain_load[bus_id,ts]
            end

            model_container.signed_flow[branch_id, ts] = branch_flow_expr
        end
    end
end

function add_flow_constraints(model_container::RSOAssessmentModel, network, TS)
    model_l = model_container.model
    for branch in Networks.get_branches(network)
        branch_id = Networks.get_id(branch)
        flow_limit_l = Networks.get_limit(branch)

        for ts in TS
            overflow_l = model_container.overflow[ts]
            signed_flow_l = model_container.signed_flow[branch_id, ts]

            name = @sprintf("c_pos_flow_limit[%s,%s]",branch_id,ts)
            @constraint(model_l, signed_flow_l <= flow_limit_l + overflow_l , base_name=name)
            name = @sprintf("c_neg_flow_limit[%s,%s]",branch_id,ts)
            @constraint(model_l, -signed_flow_l <= flow_limit_l + overflow_l , base_name=name)

        end
    end
end

function add_kkt_stationarity_constraints(model_container::RSOAssessmentModel, network, TS)
    for ts in TS
        duals_sum_l = AffExpr(0.)
        for branch_id in Networks.get_id.(Networks.get_branches(network))
            add_to_expression!(duals_sum_l,
                model_container.pos_flow_limit_duals[branch_id,ts] + model_container.neg_flow_limit_duals[branch_id,ts])  
        end
        @constraint(model_container.model, duals_sum_l == 1.)
    end

    return model_container
end

function add_kkt_complementarity_constraints(model_container, network, TS)
    for branch in Networks.get_branches(network)
        branch_id = Networks.get_id(branch)
        flow_limit_l = Networks.get_limit(branch)

        for ts in TS
            # pos_flow_limit_duals
            kkt_var = model_container.pos_flow_limit_duals[branch_id, ts]
            cstr_expr = - model_container.signed_flow[branch_id, ts] + flow_limit_l + model_container.overflow[ts]
            b_indicator = model_container.pos_flow_limit_indicators[branch_id, ts]
            ub_cstr = compute_ub(cstr_expr)
            ub_kkt = upper_bound(kkt_var)
            formulate_complementarity_constraints!(model_container.model,
                                                    kkt_var, cstr_expr, b_indicator, ub_kkt, ub_cstr)
            # neg_flow_limit_duals
            kkt_var = model_container.neg_flow_limit_duals[branch_id, ts]
            cstr_expr = model_container.signed_flow[branch_id, ts] + flow_limit_l + model_container.overflow[ts]
            b_indicator = model_container.neg_flow_limit_indicators[branch_id, ts]
            ub_cstr = compute_ub(cstr_expr)
            ub_kkt = upper_bound(kkt_var)
            formulate_complementarity_constraints!(model_container.model,
                                                    kkt_var, cstr_expr, b_indicator, ub_kkt, ub_cstr)
        end
    end

    return model_container
end

function formulate_overflow_subproblem(model_container::RSOAssessmentModel, network, TS)
    create_flow_limit_duals_and_indicators(model_container, network, TS)
    create_flow_expressions(model_container, network, TS)
    add_flow_constraints(model_container, network, TS)
    add_kkt_stationarity_constraints(model_container, network, TS)
    add_kkt_complementarity_constraints(model_container, network, TS)

    return model_container
end

function formulate_rso_assessment(network, TS, assessment_uncertainties, tso_actions, configs)
    model_container_l = RSOAssessmentModel()

    add_uncertainties_vars!(model_container_l, network, TS, assessment_uncertainties, tso_actions)
    add_market_vars!(model_container_l, network, TS, tso_actions)
    add_overflow_var!(model_container_l, TS, configs.BIG_M)

    add_use_limitables_constraints!(model_container_l, network, TS, tso_actions)
    add_cheapest_prod_constraints!(model_container_l, network, TS)
    add_eod_constraint!(model_container_l, network, TS)
    add_uncertainties_constraints!(model_container_l, network, TS, assessment_uncertainties, configs)

    formulate_overflow_subproblem(model_container_l, network, TS)

    create_objective!(model_container_l, network, TS)
    @objective(get_model(model_container_l), Max, model_container_l.full_obj)

    return model_container_l
end

function is_validated(model_container::RSOAssessmentModel)
    return value(model_container.full_obj) < 1e-09
end

function has_positive_slack(model_container::RSOAssessmentModel)::Bool
    return false # ( has_positive_value(model_container.lower.lol_model.p_loss_of_load)
           #  || has_positive_value(model_container.lower.lol_model.p_cut_prod) )
end

function get_assessment_uncertainties_lb(assessment_uncertainties, bus_or_limitable::String)
    return assessment_uncertainties[bus_or_limitable][1]
end
function get_assessment_uncertainties_ub(assessment_uncertainties, bus_or_limitable::String)
    return assessment_uncertainties[bus_or_limitable][2]
end
