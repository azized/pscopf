using .Networks

using JuMP
using Dates
using DataStructures
using Printf
using Parameters
"""
REF_SCHEDULE_TYPE_IN_TSO : Indicates which schedule to use as reference for pilotables state/levels needed
                            for sequencing constraints and TSO objective function.
"""
@with_kw mutable struct TSOBilevelConfigs
    CONSIDER_N_1_CSTRS::Bool = false
    TSO_LIMIT_PENALTY::Float64 = 1e-3
    TSO_LOL_PENALTY::Float64 = 1e5
    TSO_CAPPING_COST::Float64 = 1.
    TSO_PILOTABLE_BOUNDING_COST::Float64 = 1.
    USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST::Bool = false
    MARKET_LOL_PENALTY::Float64 = 1e5
    MARKET_CAPPING_COST::Float64 = 1.
    out_path::Union{Nothing,String} = nothing
    problem_name::String = "TSOBilevel"
    LINK_SCENARIOS_LIMIT::Bool = true
    LINK_SCENARIOS_PILOTABLE_LEVEL::Bool = false
    LINK_SCENARIOS_PILOTABLE_ON::Bool = false
    LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET::Bool = false
    big_m = 1e6
    REF_SCHEDULE_TYPE_IN_TSO::Union{Market,TSO} = Market();
end

##########################################################
#        upper problem : TSO
##########################################################

@with_kw struct TSOBilevelTSOLimitableModel <: AbstractLimitableModel
    #gen,ts,s
    p_injected = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_limit = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_is_limited = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_limit_x_is_limited = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_capping = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #ts,s
    p_global_capping = SortedDict{Tuple{DateTime,String},VariableRef}();
end

@with_kw struct TSOBilevelTSOPilotableModel <: AbstractPilotableModel
    #gen,ts,s
    p_imposition_min = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_imposition_max = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_start = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_on = SortedDict{Tuple{String,DateTime,String},VariableRef}();
end

@with_kw struct TSOBilevelTSOLoLModel <: AbstractLoLModel
    #bus,ts,s #Loss of Load
    p_loss_of_load = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #ts,s
    p_global_loss_of_load = SortedDict{Tuple{DateTime,String},VariableRef}();
end

@with_kw mutable struct TSOBilevelTSOObjectiveModel <: AbstractObjectiveModel
    limitable_cost = GenericAffExpr{Float64, VariableRef}()
    pilotable_cost = GenericAffExpr{Float64, VariableRef}()

    penalty = GenericAffExpr{Float64, VariableRef}()

    full_obj = GenericAffExpr{Float64, VariableRef}()
end

##########################################################
#        lower problem : Market
##########################################################
@with_kw struct TSOBilevelMarketLimitableModel <: AbstractLimitableModel
    #ts,s #FIXME not sure if this should be a lower or an upper variable
    p_global_injected = SortedDict{Tuple{DateTime,String},VariableRef}();
    #ts,s
    p_global_capping = SortedDict{Tuple{DateTime,String},VariableRef}();
end

@with_kw struct TSOBilevelMarketPilotableModel <: AbstractPilotableModel
    #gen,ts,s
    p_injected = SortedDict{Tuple{String,DateTime,String},VariableRef}();
end

@with_kw struct TSOBilevelMarketLoLModel <: AbstractLoLModel
    #ts,s
    p_global_loss_of_load = SortedDict{Tuple{DateTime,String},VariableRef}();
end

@with_kw mutable struct TSOBilevelMarketObjectiveModel <: AbstractObjectiveModel
    limitable_cost = GenericAffExpr{Float64, VariableRef}()
    pilotable_cost = GenericAffExpr{Float64, VariableRef}()

    penalty = GenericAffExpr{Float64, VariableRef}()

    full_obj = GenericAffExpr{Float64, VariableRef}()
end

##########################################################
#        TSOBilevel Model
##########################################################

@with_kw struct TSOBilevelTSOModelContainer <: AbstractModelContainer
    model::Model
    limitable_model::TSOBilevelTSOLimitableModel = TSOBilevelTSOLimitableModel()
    pilotable_model::TSOBilevelTSOPilotableModel = TSOBilevelTSOPilotableModel()
    lol_model::TSOBilevelTSOLoLModel = TSOBilevelTSOLoLModel()
    objective_model::TSOBilevelTSOObjectiveModel = TSOBilevelTSOObjectiveModel()
    #branch,ts,s
    flows::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
    rso_constraint::SortedDict{Tuple{String,DateTime,String},ConstraintRef} =
        SortedDict{Tuple{String,DateTime,String},ConstraintRef}()
end
@with_kw struct TSOBilevelMarketModelContainer <: AbstractModelContainer
    model::Model
    limitable_model::TSOBilevelMarketLimitableModel = TSOBilevelMarketLimitableModel()
    pilotable_model::TSOBilevelMarketPilotableModel = TSOBilevelMarketPilotableModel()
    lol_model::TSOBilevelMarketLoLModel = TSOBilevelMarketLoLModel()
    objective_model::TSOBilevelMarketObjectiveModel = TSOBilevelMarketObjectiveModel()
    #ts,s
    eod_constraint::SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef} =
        SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef}()
end
@with_kw struct TSOBilevelKKTModelContainer <: AbstractModelContainer
    #TODO create adequate struct to link each dual kkt variable, indicator variable, and constraint to each other
    # e.g.
    # eod_model = ts,s -> Struct{dual, indicator, ConstraintRef}
    # pilotable_min = id,ts,s -> Struct{dual, indicator, ConstraintRef} or Struct{dual, Nothing, ConstraintRef}
    # than simply call reformulate_kkt(cstr) when building lower problem's constraints provided that objective is created first
    model::Model
    #ts,s
    eod_duals::SortedDict{Tuple{DateTime,String},VariableRef} =
        SortedDict{Tuple{DateTime,String},VariableRef}()
    capping_duals::SortedDict{Tuple{DateTime,String},VariableRef} =
        SortedDict{Tuple{DateTime,String},VariableRef}()
    capping_indicators::SortedDict{Tuple{DateTime,String},VariableRef} =
        SortedDict{Tuple{DateTime,String},VariableRef}()
    loss_of_load_duals::SortedDict{Tuple{DateTime,String},VariableRef} =
        SortedDict{Tuple{DateTime,String},VariableRef}()
    loss_of_load_indicators::SortedDict{Tuple{DateTime,String},VariableRef} =
        SortedDict{Tuple{DateTime,String},VariableRef}()
    #pilotable_id,ts,s
    firmness_duals::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
    pmin_duals::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
    pmin_indicators::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
    pmax_duals::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
    pmax_indicators::SortedDict{Tuple{String,DateTime,String},VariableRef} =
        SortedDict{Tuple{String,DateTime,String},VariableRef}()
end
TSOBilevelModel = BilevelModelContainer{TSOBilevelTSOModelContainer, TSOBilevelMarketModelContainer, TSOBilevelKKTModelContainer}
function BilevelModelContainer{TSOBilevelTSOModelContainer,TSOBilevelMarketModelContainer,TSOBilevelKKTModelContainer}()
    bilevel_model = Model()
    upper = TSOBilevelTSOModelContainer(model=bilevel_model)
    lower = TSOBilevelMarketModelContainer(model=bilevel_model)
    kkt_model = TSOBilevelKKTModelContainer(model=bilevel_model)
    return BilevelModelContainer(bilevel_model, upper, lower, kkt_model)
end

"""
Computes the actual capping in the Limitable model

    The global capping of the upper limitable model can be lower than the actual capping
    The actual capping is given by the localised capping in the upper problem or the global capping in the lower problem
"""
function sum_capping(limitable_model::TSOBilevelTSOLimitableModel, ts,s, network::Networks.Network)
    sum_l = 0.
    for gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        sum_l += limitable_model.p_capping[gen_id,ts,s]
    end
    return sum_l
end

"""
Computes the actual loss of load (LoL) in the LoL model

    The global LoL of the upper limitable model can be lower than the actual LoL
    The actual LoL is given by the localised LoL in the upper problem or the global LoL in the lower problem
"""
function sum_lol(lol_model::TSOBilevelTSOLoLModel, ts, s, network::Networks.Network)
    sum_l = 0.
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        sum_l += lol_model.p_loss_of_load[bus_id,ts,s]
    end
    return sum_l
end

function has_positive_slack(model_container::TSOBilevelModel)::Bool
    return has_positive_value(model_container.lower.lol_model.p_global_loss_of_load) #If TSO cut => market did too
end

function get_upper_obj_expr(bilevel_model::TSOBilevelModel)
    return bilevel_model.upper.objective_model.full_obj
end
function get_lower_obj_expr(bilevel_model::TSOBilevelModel)
    return bilevel_model.lower.objective_model.full_obj
end


##########################################################
#        upper problem : TSO functions : TSOBilevelTSOModelContainer
##########################################################
function create_tso_vars!( model_container::TSOBilevelTSOModelContainer,
                            network::Networks.Network,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String}
                            )
    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    limitables_list_l = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    buses_list_l = Networks.get_buses(network)
    add_pilotables_vars!(model_container,
                        pilotables_list_l, target_timepoints, scenarios,
                        imposition_vars=true, commitment_vars=true,
                        prefix="tso")
    add_limitables_vars!(model_container,
                        target_timepoints, scenarios,
                        limitables_list_l,
                        injection_vars=true, limit_vars=true, local_capping_vars=true, global_capping_vars=true,
                        prefix="tso"
                        )
    add_lol_vars!(model_container,
                target_timepoints, scenarios, buses_list_l,
                global_lol_vars=true, local_lol_vars=true,
                prefix="tso")
end


function add_tso_constraints!(bimodel_container::TSOBilevelModel,
                            target_timepoints, scenarios, network,
                            firmness, reference_schedule, generators_initial_state,
                            uncertainties_at_ech::UncertaintiesAtEch,
                            configs::TSOBilevelConfigs)
    model = bimodel_container.model
    tso_model_container::TSOBilevelTSOModelContainer = bimodel_container.upper
    market_model_container::TSOBilevelMarketModelContainer = bimodel_container.lower

    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    unit_commitment_constraints!(model,
                    tso_model_container.pilotable_model, pilotables_list_l,  target_timepoints, scenarios,
                    firmness, reference_schedule, generators_initial_state,
                    always_link_scenarios=configs.LINK_SCENARIOS_PILOTABLE_ON,
                    prefix="tso")
    power_imposition_constraints!(model, tso_model_container.pilotable_model,
                    pilotables_list_l, target_timepoints, scenarios,
                    firmness, reference_schedule;
                    always_link_scenarios=configs.LINK_SCENARIOS_PILOTABLE_LEVEL,
                    prefix="tso")

    limitables_list_l = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    limitables_ids_l = map(lim_gen->Networks.get_id(lim_gen), limitables_list_l)
    limitable_model = tso_model_container.limitable_model
    limitable_power_constraints!(model,
                        limitable_model, limitables_list_l, target_timepoints, scenarios,
                        firmness, uncertainties_at_ech,
                        always_link_scenarios=configs.LINK_SCENARIOS_LIMIT,
                        prefix="tso")
    local_capping_constraints!(model,
                            limitable_model, limitables_list_l, target_timepoints, scenarios,
                            uncertainties_at_ech, prefix="tso")
    global_capping_constraints!(model,
                            limitable_model, limitables_list_l,
                            target_timepoints, scenarios,
                            uncertainties_at_ech,
                            prefix="tso")
    distribution_constraint!(model,
                            get_global_capping(market_model_container.limitable_model),
                            get_local_capping(tso_model_container.limitable_model),
                            limitables_ids_l, target_timepoints, scenarios,
                            cstr_prefix_name="tso_distribute_capping")
    distribution_constraint!(model,
                            market_model_container.limitable_model.p_global_injected, #TODO rename to avoid onfusion between localised injections and global ones
                            get_p_injected(tso_model_container.limitable_model),
                            limitables_ids_l, target_timepoints, scenarios,
                            cstr_prefix_name="tso_distribute_enr_injections")

    buses_list_l = Networks.get_buses(network)
    buses_ids_l = map(bus->Networks.get_id(bus), buses_list_l)
    lol_model = tso_model_container.lol_model
    local_lol_constraints!(model,
                        lol_model,
                        buses_list_l, target_timepoints, scenarios,
                        uncertainties_at_ech,
                        prefix="tso")
    global_lol_constraints!(model,
                        lol_model, buses_list_l, target_timepoints, scenarios,
                        uncertainties_at_ech,
                        prefix="tso")
    distribution_constraint!(model,
                            get_global_lol(market_model_container.lol_model),
                            get_local_lol(tso_model_container.lol_model),
                            buses_ids_l, target_timepoints, scenarios,
                            cstr_prefix_name="tso_distribute_lol")

    combinations = (configs.CONSIDER_N_1_CSTRS) ? all_combinations(network, target_timepoints, scenarios) :
                                                    all_n_combinations(network, target_timepoints, scenarios)
    rso_constraints!(bimodel_container.model,
                    tso_model_container.flows,
                    tso_model_container.rso_constraint,
                    market_model_container.pilotable_model, #pilotable injections are decided by Market
                    tso_model_container.limitable_model, #limitable injections are decided by TSO
                    tso_model_container.lol_model,
                    combinations,
                    uncertainties_at_ech, network,
                    prefix="tso")

    return bimodel_container
end

function create_tso_objectives!(model_container::TSOBilevelTSOModelContainer,
                                target_timepoints, scenarios, network,
                                preceding_market_schedule::Schedule,
                                capping_cost, loss_of_load_cost,
                                limit_penalty,
                                pilotable_bounding_cost,
                                use_prop_cost_for_bounding::Bool)
    objective_model = model_container.objective_model

    # objective_model.pilotable_cost
    add_tsobilevel_impositions_cost!(model_container,
                                    target_timepoints, scenarios, network,
                                    preceding_market_schedule,
                                    pilotable_bounding_cost,
                                    use_prop_cost_for_bounding)

    # limitable_cost : capping (fr. ecretement)
    objective_model.limitable_cost += coeffxsum(model_container.limitable_model.p_global_capping, capping_cost)

    # cost for cutting consumption (lol) and avoid limiting for no reason
    objective_model.penalty += coeffxsum(model_container.lol_model.p_global_loss_of_load, loss_of_load_cost)
    objective_model.penalty += coeffxsum(model_container.limitable_model.b_is_limited, limit_penalty)

    objective_model.full_obj = ( objective_model.pilotable_cost +
                                objective_model.limitable_cost +
                                objective_model.penalty )
    @objective(model_container.model, Min, objective_model.full_obj)
    return model_container
end

function add_tsobilevel_impositions_cost!(model_container::TSOBilevelTSOModelContainer,
                                        target_timepoints, scenarios, network,
                                        preceding_market_schedule,
                                        pilotable_bounding_cost,
                                        use_prop_cost_for_bounding::Bool)
    tso_pilotable_model = model_container.pilotable_model
    objective_expr = model_container.objective_model.pilotable_cost

    for gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        gen_id = Networks.get_id(gen)
        if use_prop_cost_for_bounding
            cost = Networks.get_prop_cost(gen)
        else
            cost = pilotable_bounding_cost
        end

        for ts in target_timepoints
            for s in scenarios
                commitment = get_commitment_value(preceding_market_schedule, gen_id, ts, s)
                if  !Networks.needs_commitment(gen) || (!ismissing(commitment) && (commitment==ON))
                    add_tsobilevel_started_impositions_cost!(objective_expr, gen,
                                                            tso_pilotable_model.p_imposition_min[gen_id, ts, s],
                                                            tso_pilotable_model.p_imposition_max[gen_id, ts, s],
                                                            cost)
                else
                    add_tsobilevel_non_started_impositions_cost!(objective_expr, model_container,
                                                                 gen, ts, s, cost)
                end
            end
        end
    end

    return model_container
end
function add_tsobilevel_started_impositions_cost!(objective_expr::AffExpr,
                                                gen::Generator,
                                                pmin_var, pmax_var,
                                                pilotable_bounding_cost
                                                )
    p_max = Networks.get_p_max(gen)
    p_min = Networks.get_p_min(gen)

    add_to_expression!(objective_expr, pilotable_bounding_cost * (p_max - pmax_var))
    add_to_expression!(objective_expr, pilotable_bounding_cost * (pmin_var - p_min))

    return objective_expr
end
function add_tsobilevel_non_started_impositions_cost!(objective_expr::AffExpr,
                                                      tso_model_container::TSOBilevelTSOModelContainer,
                                                    gen, ts, s,
                                                    pilotable_bounding_cost
                                                    )
    #need a second expression for units that do not need a commitment (p_min=0 and no b_on)
    @assert Networks.needs_commitment(gen)

    tso_pilotable_model = tso_model_container.pilotable_model
    gen_id = Networks.get_id(gen)

    p_max = Networks.get_p_max(gen)
    # p_min = Networks.get_p_min(gen)

    #need to cost reducing pmax otherwise TSO may always limit pmax when starting a unit
    pmax_var = tso_pilotable_model.p_imposition_max[gen_id, ts, s]
    b_on_var = tso_pilotable_model.b_on[gen_id, ts, s]
    add_to_expression!(objective_expr, pilotable_bounding_cost * (p_max*b_on_var - pmax_var) )

    pmin_var = tso_pilotable_model.p_imposition_min[gen_id, ts, s]
    add_to_expression!(objective_expr, pilotable_bounding_cost * pmin_var)
    # add_to_expression!(objective_expr, pilotable_bounding_cost * p_min * b_on_var)

    return objective_expr
end

##########################################################
#        lower problem : Market functtons  : TSOBilevelMarketModelContainer
##########################################################

function create_market_vars!(model_container::TSOBilevelMarketModelContainer,
                            network::Networks.Network,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String}
                            )

    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    add_pilotables_vars!(model_container,
                        pilotables_list_l, target_timepoints, scenarios,
                        injection_vars=true,
                        prefix="market_"
                        )
    add_limitables_vars!(model_container,
                        target_timepoints, scenarios,
                        global_capping_vars=true,
                        global_injection_vars=true,
                        prefix="market_"
                        )
    add_lol_vars!(model_container,
                target_timepoints, scenarios,
                global_lol_vars=true,
                prefix="market_")
end


function add_firmness_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                            pilotables_list, target_timepoints, scenarios,
                            firmness,
                            link_scenarios_pilotable_level_market;
                            prefix::String="")
    for gen_pilotable in pilotables_list
        gen_id = Networks.get_id(gen_pilotable)
        for ts in target_timepoints
            decision_firmness_l = get_power_level_firmness(firmness, gen_id, ts)
            for s in scenarios[2:end]
                if requires_linking(decision_firmness_l, link_scenarios_pilotable_level_market)
                    #create duals relative to firmness constraints of injected pilotable power in market
                    name = @sprintf("c_firmness[%s,%s,%s]",gen_id,ts,s)
                    add_dual!(kkt_model_container.model, kkt_model_container.firmness_duals,
                                (gen_id,ts,s), name, false)
                end
            end
        end
    end
    return kkt_model_container
end

function add_eod_constraints_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                            target_timepoints, scenarios)
    for ts in target_timepoints
        for s in scenarios
            #create duals relative to EOD constraint
            name = @sprintf("c_eod[%s,%s]",ts,s)
            add_dual!(kkt_model_container.model, kkt_model_container.eod_duals, (ts,s), name, false)
        end
    end
end


function add_min_global_capping_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                                    target_timepoints, scenarios)
    for ts in target_timepoints
        for s in scenarios
            #create duals and indicators relative to TSO min capping constraint
            name = @sprintf("c_min_e[%s,%s]",ts,s)
            add_dual_and_indicator!(kkt_model_container.model,
                                    kkt_model_container.capping_duals, kkt_model_container.capping_indicators, (ts,s),
                                    name, true)
        end
    end
end


function add_min_global_lol_duals!(market_model_container::TSOBilevelMarketModelContainer,
                                kkt_model_container::TSOBilevelKKTModelContainer,
                                target_timepoints, scenarios)
    for ts in target_timepoints
        for s in scenarios
            #create duals and indicators relative to TSO min LoL constraint
            name = @sprintf("c_min_lol[%s,%s]",ts,s)
            add_dual_and_indicator!(kkt_model_container.model,
                                    kkt_model_container.loss_of_load_duals, kkt_model_container.loss_of_load_indicators, (ts,s),
                                    name, true)
        end
    end
    return market_model_container
end
function add_impositions_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                            pilotables_list, target_timepoints, scenarios)
    for ts in target_timepoints
        for s in scenarios
            for gen in pilotables_list
                gen_id = Networks.get_id(gen)
                add_pilotable_pmin_duals!(kkt_model_container,
                                        gen_id, ts, s)
                add_pilotable_pmax_duals!(kkt_model_container,
                                        gen_id, ts, s)
            end
        end
    end

    return kkt_model_container
end
function add_pilotable_pmin_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                                gen_id, ts, s)
    #create duals and indicators relative to tso pmin constraint
    name = @sprintf("c_tso_pmin[%s,%s,%s]",gen_id,ts,s)
    add_dual_and_indicator!(kkt_model_container.model,
                            kkt_model_container.pmin_duals, kkt_model_container.pmin_indicators, (gen_id,ts,s),
                            name, true)

    return kkt_model_container
end
function add_pilotable_pmax_duals!(kkt_model_container::TSOBilevelKKTModelContainer,
                                gen_id, ts, s)
    #create duals and indicators relative to tso pmax constraint
    name = @sprintf("c_tso_pmax[%s,%s,%s]",gen_id,ts,s)
    add_dual_and_indicator!(kkt_model_container.model,
                            kkt_model_container.pmax_duals, kkt_model_container.pmax_indicators, (gen_id,ts,s),
                            name, true)

    return kkt_model_container
end


function add_market_constraints!(bimodel_container::TSOBilevelModel,
                            target_timepoints, scenarios, network,
                            firmness, #reference_schedule, generators_initial_state,
                            uncertainties_at_ech::UncertaintiesAtEch,
                            configs::TSOBilevelConfigs)
    model = bimodel_container.model
    tso_model_container::TSOBilevelTSOModelContainer = bimodel_container.upper
    market_model_container::TSOBilevelMarketModelContainer = bimodel_container.lower
    kkt_model_container::TSOBilevelKKTModelContainer = bimodel_container.kkt_model

    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    pilotable_power_constraints!(model,
                                market_model_container.pilotable_model,
                                pilotables_list_l, target_timepoints, scenarios,
                                firmness,
                                missing,
                                ignore_sequencing_cstrs=true, #Sequencing constraints are imposed in TSO upper problem
                                always_link_scenarios=configs.LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET,
                                prefix="market_")
    add_firmness_duals!(kkt_model_container,
                        pilotables_list_l, target_timepoints, scenarios, firmness,
                        configs.LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET,
                        prefix="market_")
    respect_impositions_constraints!(model,
                                market_model_container.pilotable_model, #Injections are in the lower problem
                                tso_model_container.pilotable_model, #Impositions are in the upper problem
                                pilotables_list_l,  target_timepoints, scenarios,
                                prefix="market_")
    add_impositions_duals!(kkt_model_container,
                        pilotables_list_l, target_timepoints, scenarios)

    limitables_list_l = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    global_capping_constraints!(model,
                            market_model_container.limitable_model,
                            limitables_list_l, target_timepoints, scenarios,
                            uncertainties_at_ech,
                            min_cap=get_global_capping(tso_model_container.limitable_model)
                            )
    add_min_global_capping_duals!(kkt_model_container,
                            target_timepoints, scenarios)
    global_injected_constraints!(model,
                            market_model_container.limitable_model,
                            limitables_list_l, target_timepoints, scenarios,
                            uncertainties_at_ech, prefix="market_")

    buses_list_l = Networks.get_buses(network)
    lol_model = market_model_container.lol_model
    global_lol_constraints!(model,
                        lol_model, buses_list_l, target_timepoints, scenarios,
                        uncertainties_at_ech,
                        min_lol=get_global_lol(tso_model_container.lol_model),
                        prefix="market")
    add_min_global_lol_duals!(market_model_container, kkt_model_container,
                        target_timepoints, scenarios)

    eod_constraints!(market_model_container.model, market_model_container.eod_constraint,
                    market_model_container.pilotable_model,
                    tso_model_container.limitable_model,
                    market_model_container.lol_model,
                    target_timepoints, scenarios,
                    uncertainties_at_ech, network
                    )
    add_eod_constraints_duals!(kkt_model_container, target_timepoints, scenarios)

    return bimodel_container
end

function create_market_objectives!(model_container::TSOBilevelMarketModelContainer,
                                network,
                                capping_cost, loss_of_load_cost)
    objective_model = model_container.objective_model

    # model_container.pilotable_cost
    add_prop_cost!(model_container.objective_model.pilotable_cost,
                            model_container.pilotable_model.p_injected, network)

    # limitable_cost : capping (fr. ecretement)
    objective_model.limitable_cost += coeffxsum(model_container.limitable_model.p_global_capping, capping_cost)

    # cost for cutting load/consumption
    objective_model.penalty += coeffxsum(model_container.lol_model.p_global_loss_of_load, loss_of_load_cost)

    objective_model.full_obj = ( objective_model.pilotable_cost +
                                objective_model.limitable_cost +
                                objective_model.penalty )

    return model_container
end

##########################################################
#        kkt reformulation : TSOBilevelKKTModelContainer
##########################################################

function add_kkt_stationarity_constraints!(kkt_model::TSOBilevelKKTModelContainer,
                                            target_timepoints, scenarios, network,
                                            capping_cost, loss_of_load_cost)
    #FIXME can be generic by iterating on lower variables to construct each stationarity constraint
    # iterate on the objective and lower constraints to extract their coefficients, but need to link each cnstraint to its dual var
    add_capping_stationarity_constraints!(kkt_model, target_timepoints, scenarios, capping_cost)
    add_loss_of_load_stationarity_constraints!(kkt_model, target_timepoints, scenarios, loss_of_load_cost)
    add_pilotable_stationarity_constraints!(kkt_model, target_timepoints, scenarios, network)
end

function add_loss_of_load_stationarity_constraints!(kkt_model::TSOBilevelKKTModelContainer,
                                                target_timepoints, scenarios, loss_of_load_cost)
    for ts in target_timepoints
        for s in scenarios
            # @assert ( capping_cost ≈ coefficient(market_model.objective_model.full_obj,
            #                                                 market_model.lol_model.p_global_loss_of_load[ts,s]) )
            name = @sprintf("c_stationarity_lol[%s,%s]",ts,s)
            @constraint(kkt_model.model,
                        loss_of_load_cost + kkt_model.eod_duals[ts,s] - kkt_model.loss_of_load_duals[ts,s] == 0,
                        base_name=name)
        end
    end
end

function add_capping_stationarity_constraints!(kkt_model::TSOBilevelKKTModelContainer,
                                                target_timepoints, scenarios, capping_cost)
    for ts in target_timepoints
        for s in scenarios
            # @assert ( capping_cost ≈ coefficient(market_model.objective_model.full_obj,
            #                                                 market_model.limitable_model.p_global_capping[ts,s]) )
            name = @sprintf("c_stationarity_e[%s,%s]",ts,s)
            @constraint(kkt_model.model,
                        capping_cost - kkt_model.eod_duals[ts,s] - kkt_model.capping_duals[ts,s] == 0,
                        base_name=name)
        end
    end
end

function firmness_duals_sationarity_expr(kkt_model,
                                        gen_id, ts, scenario,
                                        scenarios)::AffExpr
    if length(scenarios) == 1
        result_l = 0.
    elseif scenario == scenarios[1]
        sum_l = sum( get(kkt_model.firmness_duals, (gen_id,ts,s), 0.) for s in scenarios[2:end])
        result_l = -sum_l
    else
        result_l = get(kkt_model.firmness_duals, (gen_id,ts,scenario), 0.)
    end
    return result_l
end
function add_pilotable_stationarity_constraints!(kkt_model::TSOBilevelKKTModelContainer,
                                                target_timepoints, scenarios, network)
    for pilotable_gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        gen_id = Networks.get_id(pilotable_gen)
        gen_prop_cost = Networks.get_prop_cost(pilotable_gen)
        for ts in target_timepoints
            for s in scenarios
                # @assert ( pilotable_bounding_cost ≈ coefficient(market_model.objective_model.full_obj,
                #                                             market_model.pilotable_model.p_injected[gen_id,ts,s]) )
                name = @sprintf("c_stationarity_pilotable_p[%s,%s,%s]",gen_id,ts,s)
                firmness_sationarity_expr = firmness_duals_sationarity_expr(kkt_model,
                                                                                gen_id, ts, s,
                                                                                scenarios)
                @constraint(kkt_model.model,
                            0 == gen_prop_cost
                                    + kkt_model.eod_duals[ts,s]
                                    - kkt_model.pmin_duals[gen_id,ts,s]
                                    + kkt_model.pmax_duals[gen_id,ts,s]
                                    + firmness_sationarity_expr,
                            base_name=name)
            end
        end
    end
end

function add_kkt_complementarity_constraints!(model_container::TSOBilevelModel,
                                            big_m, target_timepoints, scenarios, network)
    #FIXME can be done iteratively and generically if we loop on constraint expressions and know their corresponding kkt vars
    add_emin_complementarity_constraints!(model_container, big_m, target_timepoints, scenarios)
    add_lolmin_complementarity_constraints!(model_container, big_m, target_timepoints, scenarios)
    add_pmin_complementarity_constraints!(model_container, big_m, target_timepoints, scenarios, network)
    add_pmax_complementarity_constraints!(model_container, big_m, target_timepoints, scenarios, network)
end

"""
    complementarity constraint linked to @ref PSCOPF.global_capping_constraints!
"""
function add_emin_complementarity_constraints!(model_container::TSOBilevelModel,
                                            big_m, target_timepoints, scenarios)
    tso_limitable_model = model_container.upper.limitable_model
    market_limitable_model = model_container.lower.limitable_model
    kkt_model = model_container.kkt_model

    for ts in target_timepoints
        for s in scenarios
            kkt_var = kkt_model.capping_duals[ts,s]
            cstr_expr = market_limitable_model.p_global_capping[ts,s] - tso_limitable_model.p_global_capping[ts,s]
            b_indicator = kkt_model.capping_indicators[ts,s]
            ub_cstr = compute_ub(cstr_expr, big_m)
            formulate_complementarity_constraints!(kkt_model.model, kkt_var, cstr_expr, b_indicator, big_m, ub_cstr)
        end
    end
    return model_container
end

"""
    complementarity constraint linked to @ref PSCOPF.global_lol_constraints!
"""
function add_lolmin_complementarity_constraints!(model_container::TSOBilevelModel,
                                                big_m, target_timepoints, scenarios)
    tso_lol_model = model_container.upper.lol_model
    market_lol_model = model_container.lower.lol_model
    kkt_model = model_container.kkt_model

    for ts in target_timepoints
        for s in scenarios
            kkt_var = kkt_model.loss_of_load_duals[ts,s]
            cstr_expr = market_lol_model.p_global_loss_of_load[ts,s] - tso_lol_model.p_global_loss_of_load[ts,s]
            b_indicator = kkt_model.loss_of_load_indicators[ts,s]
            ub_cstr = compute_ub(cstr_expr, big_m)
            formulate_complementarity_constraints!(kkt_model.model, kkt_var, cstr_expr, b_indicator, big_m, ub_cstr)
        end
    end
    return model_container
end

"""
    complementarity constraint linked to the p_min part of @ref PSCOPF.respect_impositions_constraints!
"""
function add_pmin_complementarity_constraints!(model_container::TSOBilevelModel,
                                                big_m, target_timepoints, scenarios, network)
    tso_pilotable_model = model_container.upper.pilotable_model
    market_pilotable_model = model_container.lower.pilotable_model
    kkt_model = model_container.kkt_model

    for ts in target_timepoints
        for s in scenarios
            for pilotable_gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
                gen_id = Networks.get_id(pilotable_gen)

                kkt_var = kkt_model.pmin_duals[gen_id,ts,s]
                cstr_expr = market_pilotable_model.p_injected[gen_id,ts,s] - tso_pilotable_model.p_imposition_min[gen_id,ts,s]
                b_indicator = kkt_model.pmin_indicators[gen_id,ts,s]
                ub_cstr = compute_ub(cstr_expr, big_m) #or get_p_max(pilotable_gen)
                formulate_complementarity_constraints!(kkt_model.model, kkt_var, cstr_expr, b_indicator, big_m, ub_cstr)
            end
        end
    end
    return model_container
end

"""
    complementarity constraint linked to the p_max part of @ref PSCOPF.respect_impositions_constraints!
"""
function add_pmax_complementarity_constraints!(model_container::TSOBilevelModel,
                                                big_m, target_timepoints, scenarios, network)
    tso_pilotable_model = model_container.upper.pilotable_model
    market_pilotable_model = model_container.lower.pilotable_model
    kkt_model = model_container.kkt_model

    for ts in target_timepoints
        for s in scenarios
            for pilotable_gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
                gen_id = Networks.get_id(pilotable_gen)

                kkt_var = kkt_model.pmax_duals[gen_id,ts,s]
                cstr_expr = tso_pilotable_model.p_imposition_max[gen_id,ts,s] - market_pilotable_model.p_injected[gen_id,ts,s]
                b_indicator = kkt_model.pmax_indicators[gen_id,ts,s]
                ub_cstr = compute_ub(cstr_expr, big_m) #or get_p_max(pilotable_gen)
                formulate_complementarity_constraints!(kkt_model.model, kkt_var, cstr_expr, b_indicator, big_m, ub_cstr)
            end
        end
    end
    return model_container
end


#############################
# Utils
#############################
function coeffxsum(vars_dict::AbstractDict{T,V}, coeff::Float64
                )::GenericAffExpr{Float64, VariableRef} where T <: Tuple where V <: AbstractVariableRef
    terms = [var_l=>coeff for (_, var_l) in vars_dict]
    result = GenericAffExpr{Float64, VariableRef}(0., terms)
    return result
end

##########################################################
#        TSOBilevel
##########################################################
function tso_bilevel(network::Networks.Network,
                    target_timepoints::Vector{Dates.DateTime},
                    generators_initial_state::SortedDict{String,GeneratorState},
                    scenarios::Vector{String},
                    uncertainties_at_ech::UncertaintiesAtEch,
                    firmness::Firmness,
                    preceding_market_schedule::Schedule,
                    preceding_tso_schedule::Schedule,
                    # gratis_starts::Set{Tuple{String,Dates.DateTime}},
                    configs::TSOBilevelConfigs
                    )

    @assert(configs.big_m >= configs.MARKET_LOL_PENALTY)
    @assert(configs.big_m >= configs.MARKET_CAPPING_COST)
    @assert(all( configs.big_m >= Networks.get_prop_cost(gen)
                for gen in Networks.get_generators_of_type(network, Networks.PILOTABLE) ))

    if is_market(configs.REF_SCHEDULE_TYPE_IN_TSO)
        reference_schedule = preceding_market_schedule
    elseif is_tso(configs.REF_SCHEDULE_TYPE_IN_TSO)
        reference_schedule = preceding_tso_schedule
    else
        throw( error("Invalid REF_SCHEDULE_TYPE_IN_TSO config.") )
    end

    bimodel_container_l = TSOBilevelModel()

    create_tso_vars!(bimodel_container_l.upper,
                    network, target_timepoints, scenarios)
    create_market_vars!(bimodel_container_l.lower,
                        network, target_timepoints, scenarios)

    #this is the expression no objective is added to the jump model
    create_market_objectives!(bimodel_container_l.lower, network,
                            configs.MARKET_CAPPING_COST, configs.MARKET_LOL_PENALTY)

    #constraints may use upper and lower vars at the same time
    add_tso_constraints!(bimodel_container_l, target_timepoints, scenarios, network,
                        firmness, reference_schedule, generators_initial_state,
                        uncertainties_at_ech, configs)
    add_market_constraints!(bimodel_container_l, target_timepoints, scenarios, network, firmness, uncertainties_at_ech, configs)
    #kkt stationarity
    add_kkt_stationarity_constraints!(bimodel_container_l.kkt_model,
                                    target_timepoints, scenarios, network,
                                    configs.MARKET_CAPPING_COST, configs.MARKET_LOL_PENALTY)
    #kkt complementarity
    add_kkt_complementarity_constraints!(bimodel_container_l, configs.big_m, target_timepoints, scenarios, network)

    create_tso_objectives!(bimodel_container_l.upper,
                        target_timepoints, scenarios, network,
                        reference_schedule, #reference to see which units are currently on
                        configs.TSO_CAPPING_COST, configs.TSO_LOL_PENALTY,
                        configs.TSO_LIMIT_PENALTY,
                        configs.TSO_PILOTABLE_BOUNDING_COST, configs.USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST)

    launch_solve!(bimodel_container_l, configs)

    return bimodel_container_l
end

function launch_solve!(bimodel_container::TSOBilevelModel, configs::TSOBilevelConfigs)
    # first_iter = true
    # while(!isempty(to_add) || first_iter)
    #     first_iter = false
    #     solve!(bimodel_container, configs.problem_name, configs.out_path)

    #     to_add = verify_rso()
    #     add_constraints(to_add)
    # end

    solve!(bimodel_container, configs.problem_name, configs.out_path)
    @info("Lower Objective Value : $(value(bimodel_container.lower.objective_model.full_obj))")
end
