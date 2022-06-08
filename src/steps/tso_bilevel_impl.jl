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
    p_injected = SortedDict{Tuple{DateTime,String},VariableRef}();
    #ts,s
    p_capping = SortedDict{Tuple{DateTime,String},VariableRef}();
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

function sum_capping(limitable_model::TSOBilevelTSOLimitableModel, ts,s, network::Networks.Network)
    sum_l = 0.
    for gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        sum_l += limitable_model.p_capping[gen_id,ts,s]
    end
    return sum_l
end

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
                            generators_initial_state::SortedDict{String,GeneratorState},
                            scenarios::Vector{String},
                            uncertainties_at_ech::UncertaintiesAtEch,
                            firmness::Firmness,
                            reference_schedule::Schedule,
                            configs::TSOBilevelConfigs)
    add_limitables!(model_container,
                            network, target_timepoints, scenarios,
                            uncertainties_at_ech, firmness,
                            configs.LINK_SCENARIOS_LIMIT)
    add_pilotables!(model_container,
                            network, target_timepoints, scenarios,
                            generators_initial_state,
                            firmness,
                            reference_schedule,
                            configs.LINK_SCENARIOS_PILOTABLE_ON, configs.LINK_SCENARIOS_PILOTABLE_LEVEL)
    add_slacks!(model_container, network, target_timepoints, scenarios, uncertainties_at_ech)
end

function add_limitable!(limitable_model::TSOBilevelTSOLimitableModel, model::AbstractModel,
                            generator::Networks.Generator,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String},
                            inject_uncertainties::InjectionUncertainties,
                            power_level_firmness::SortedDict{Dates.DateTime, DecisionFirmness},#by ts
                            always_link_scenarios_limit
                        )
    gen_id = Networks.get_id(generator)
    gen_pmax = Networks.get_p_max(generator)
    for ts in target_timepoints
        for s in scenarios
            p_uncert = get_uncertainties(inject_uncertainties, ts, s)
            p_enr = min(gen_pmax, p_uncert)
            p_inj_var = add_p_injected!(limitable_model, model, gen_id, ts, s, p_enr, false)
            name =  @sprintf("P_capping[%s,%s,%s]", gen_id, ts, s)
            limitable_model.p_capping[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=p_enr)
            name =  @sprintf("c_define_e[%s,%s,%s]", gen_id, ts, s)
            @constraint(model, limitable_model.p_capping[gen_id, ts, s] == p_uncert - p_inj_var, base_name=name)
        end
        add_p_limit!(limitable_model, model, gen_id, ts, scenarios, gen_pmax,
                inject_uncertainties,
                power_level_firmness[ts],
                always_link_scenarios_limit)
    end
    return limitable_model
end
function add_limitables!(model_container::TSOBilevelTSOModelContainer,
                            network::Networks.Network,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String},
                            uncertainties_at_ech::UncertaintiesAtEch,
                            firmness::Firmness,
                            always_link_scenarios_limit=false
                            )
    model = model_container.model
    limitable_model = model_container.limitable_model
    for generator in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(generator)
        inject_uncertainties = get_uncertainties(uncertainties_at_ech, gen_id)
        add_limitable!(limitable_model, model,
                        generator, target_timepoints, scenarios,
                        inject_uncertainties,
                        get_power_level_firmness(firmness, gen_id), always_link_scenarios_limit)
    end

    for ts in target_timepoints
        for s in scenarios
            enr_max = compute_prod(uncertainties_at_ech, network, ts, s)
            name =  @sprintf("P_capping_min[%s,%s]", ts, s)
            limitable_model.p_global_capping[ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=enr_max)
        end
    end

    return model_container
end

function add_injection_bounds_sequencing_constraints!(pilotable_model, model,
                                                gen_id, target_timepoints, scenarios,
                                                gen_power_firmness,
                                                generator_reference_schedule::GeneratorSchedule,
                                                )
    # Sequencing constraints need to be applied on TSO to force them in impositions
    #to assure they are transmitted to the next market step

    for ts in target_timepoints
        if gen_power_firmness[ts] in [TO_DECIDE, DECIDED]

            if gen_power_firmness[ts] == DECIDED
            # We're past DP => p_imposition_min(s) == p_imposition_max(s) == already_decided_level \forall s
                imposed_value = safeget_prod_value(generator_reference_schedule, ts)
                freeze_vars!(model, pilotable_model.p_imposition_min, gen_id, ts, scenarios, imposed_value,
                            name="c_decided_level_p_tso_min")
                freeze_vars!(model, pilotable_model.p_imposition_max, gen_id, ts, scenarios, imposed_value,
                            name="c_decided_level_p_tso_max")
                @debug @sprintf("imposed prod level [%s,%s] : %s", gen_id, ts, imposed_value)
            end
        end
    end
end
function create_injection_bounds_vars!(pilotable_model, model,
                                    generator, target_timepoints, scenarios,
                                    gen_power_firmness::SortedDict{Dates.DateTime, DecisionFirmness},
                                    generator_reference_schedule::GeneratorSchedule,
                                    always_link_scenarios)
    gen_id = Networks.get_id(generator)
    p_max = Networks.get_p_max(generator)
    for ts in target_timepoints
        for s in scenarios
            name =  @sprintf("P_tso_min[%s,%s,%s]", gen_id, ts, s)
            pilotable_model.p_imposition_min[gen_id, ts, s] = @variable(model, base_name=name,
                                                                lower_bound=0., upper_bound=p_max)
            name =  @sprintf("P_tso_max[%s,%s,%s]", gen_id, ts, s)
            pilotable_model.p_imposition_max[gen_id, ts, s] = @variable(model, base_name=name,
                                                                lower_bound=0., upper_bound=p_max)
        end
    end

    add_scenarios_linking_constraints!(model, generator,
                                        pilotable_model.p_imposition_min,
                                        target_timepoints, scenarios,
                                        gen_power_firmness,
                                        always_link_scenarios
                                        )
    add_scenarios_linking_constraints!(model, generator,
                                        pilotable_model.p_imposition_max,
                                        target_timepoints, scenarios,
                                        gen_power_firmness,
                                        always_link_scenarios
                                        )

    add_injection_bounds_sequencing_constraints!(pilotable_model, model,
                                                gen_id, target_timepoints, scenarios,
                                                gen_power_firmness,
                                                generator_reference_schedule
                                                )

end
function create_commitment_vars!(pilotable_model::TSOBilevelTSOPilotableModel, model::Model,
                                generator, target_timepoints, scenarios,
                                generator_initial_state,
                                gen_commitment_firmness,
                                generator_reference_schedule,
                                always_link_scenarios)
    b_on_vars = pilotable_model.b_on
    b_start_vars = pilotable_model.b_start

    gen_id = Networks.get_id(generator)

    for s in scenarios
        for ts in target_timepoints
            name =  @sprintf("B_on[%s,%s,%s]", gen_id, ts, s)
            b_on_vars[gen_id, ts, s] = @variable(model, base_name=name, binary=true)
            name =  @sprintf("B_start[%s,%s,%s]", gen_id, ts, s)
            b_start_vars[gen_id, ts, s] = @variable(model, base_name=name, binary=true)
        end
    end

    #link b_on and b_start
    add_commitment_constraints!(model,
                                b_on_vars, b_start_vars,
                                gen_id, target_timepoints, scenarios, generator_initial_state)
    #linking b_on scenarios => linking b_start
    add_scenarios_linking_constraints!(model,
                                generator, b_on_vars,
                                target_timepoints, scenarios,
                                gen_commitment_firmness, always_link_scenarios
                                )
    #DMO related constraints
    add_commitment_sequencing_constraints!(model, generator,
                                        b_on_vars, b_start_vars,
                                        target_timepoints, scenarios,
                                        gen_commitment_firmness,
                                        generator_reference_schedule
                                        )
end
function add_pilotable!(pilotable_model::TSOBilevelTSOPilotableModel, model::AbstractModel,
                            generator::Networks.Generator,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String},
                            generator_initial_state::GeneratorState,
                            commitment_firmness::Union{Missing,SortedDict{Dates.DateTime, DecisionFirmness}}, #by ts #or Missing
                            power_level_firmness::SortedDict{Dates.DateTime, DecisionFirmness}, #by ts
                            generator_reference_schedule::GeneratorSchedule,
                            always_link_commitment::Bool,
                            always_link_levels::Bool)
    create_injection_bounds_vars!(pilotable_model, model, generator, target_timepoints, scenarios,
                                power_level_firmness,
                                generator_reference_schedule,
                                always_link_levels)
    if Networks.needs_commitment(generator)
        create_commitment_vars!(pilotable_model, model,
                                generator, target_timepoints, scenarios,
                                generator_initial_state,
                                commitment_firmness,
                                generator_reference_schedule,
                                always_link_commitment)
    end
end
function add_pilotables!(model_container::TSOBilevelTSOModelContainer,
                                network::Networks.Network,
                                target_timepoints::Vector{Dates.DateTime},
                                scenarios::Vector{String},
                                generators_initial_state::SortedDict{String,GeneratorState},
                                firmness::Firmness,
                                preceding_tso_schedule::Schedule,
                                always_link_commitment::Bool=false,
                                always_link_levels::Bool=false
                                )
    model = model_container.model
    pilotable_model = model_container.pilotable_model
    for pilotable_gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        gen_id = Networks.get_id(pilotable_gen)
        gen_initial_state = get_initial_state(generators_initial_state, pilotable_gen)
        commitment_firmness = get_commitment_firmness(firmness, gen_id)
        power_level_firmness = get_power_level_firmness(firmness, gen_id)
        generator_reference_schedule = get_sub_schedule(preceding_tso_schedule, gen_id)
        add_pilotable!(pilotable_model, model,
                        pilotable_gen, target_timepoints, scenarios,
                        gen_initial_state, commitment_firmness, power_level_firmness,
                        generator_reference_schedule,
                        always_link_commitment, always_link_levels)
    end
end

function add_slacks!(model_container::TSOBilevelTSOModelContainer,
                            network::Networks.Network,
                            target_timepoints::Vector{Dates.DateTime},
                            scenarios::Vector{String},
                            uncertainties_at_ech::UncertaintiesAtEch)
    #TODO: for now same as add_slacks!(::TSOModel,...)
    model = model_container.model
    lol_model = model_container.lol_model
    p_loss_of_load = lol_model.p_loss_of_load
    p_global_loss_of_load = lol_model.p_global_loss_of_load

    for ts in target_timepoints
        for s in scenarios
            conso_max = compute_load(uncertainties_at_ech, network, ts, s)
            name =  @sprintf("P_loss_of_load_min[%s,%s]", ts, s)
            p_global_loss_of_load[ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=conso_max)
        end
    end

    buses = Networks.get_buses(network)
    add_loss_of_load_by_bus!(model, p_loss_of_load,
                        buses, target_timepoints, scenarios, uncertainties_at_ech)

    return lol_model
end


function add_loss_of_load_distribution_constraint!(tso_model_container::TSOBilevelTSOModelContainer,
                                            market_lol_model::TSOBilevelMarketLoLModel,
                                            target_timepoints, scenarios, network)
    tso_model = tso_model_container.model
    tso_lol_model = tso_model_container.lol_model
    buses_ids = Networks.get_id.(Networks.get_buses(network))
    for ts in target_timepoints
        for s in scenarios
            name = @sprintf("c_dist_LOL[%s,%s]",ts,s)
            vars_sum = sum(tso_lol_model.p_loss_of_load[bus_id, ts, s]
                            for bus_id in buses_ids)
            @constraint(tso_model, market_lol_model.p_global_loss_of_load[ts, s] == vars_sum, base_name=name)
        end
    end
    return tso_model_container
end

function add_enr_distribution_constraint!(tso_model_container::TSOBilevelTSOModelContainer,
                                        market_limitable_model::TSOBilevelMarketLimitableModel,
                                        target_timepoints, scenarios, network)
    tso_model = tso_model_container.model
    limitable_model = tso_model_container.limitable_model
    limitables_ids = Networks.get_id.(Networks.get_generators_of_type(network, Networks.LIMITABLE))
    if !isempty(limitables_ids)
        for ts in target_timepoints
            for s in scenarios
                name = @sprintf("c_dist_penr[%s,%s]",ts,s)
                vars_sum = sum(limitable_model.p_injected[gen_id, ts, s]
                                for gen_id in limitables_ids)
                @constraint(tso_model, market_limitable_model.p_injected[ts, s] == vars_sum, base_name=name)
            end
        end
    end
    return tso_model_container
end

function add_capping_distribution_constraint!(tso_model_container::TSOBilevelTSOModelContainer,
                                            market_limitable_model::TSOBilevelMarketLimitableModel,
                                            target_timepoints, scenarios, network)
    tso_model = tso_model_container.model
    tso_limitable_model = tso_model_container.limitable_model
    limitables_ids = Networks.get_id.(Networks.get_generators_of_type(network, Networks.LIMITABLE))
    if !isempty(limitables_ids)
        for ts in target_timepoints
            for s in scenarios
                name = @sprintf("c_dist_e[%s,%s]",ts,s)
                vars_sum = sum(tso_limitable_model.p_capping[gen_id, ts, s]
                                for gen_id in limitables_ids)
                @constraint(tso_model, market_limitable_model.p_capping[ts, s] == vars_sum, base_name=name)
            end
        end
    end
    return tso_model_container
end
function add_injection_commitment_constraints!(tso_model_container::TSOBilevelTSOModelContainer,
                                            market_pilotable_model::TSOBilevelMarketPilotableModel,
                                            target_timepoints, scenarios, network)
    tso_model = tso_model_container.model
    b_on_vars = tso_model_container.pilotable_model.b_on
    p_imposition_min = tso_model_container.pilotable_model.p_imposition_min
    p_imposition_max = tso_model_container.pilotable_model.p_imposition_max

    p_market_inj = market_pilotable_model.p_injected

    for generator in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        if Networks.needs_commitment(generator)
            gen_id = Networks.get_id(generator)
            p_max = Networks.get_p_max(generator)
            p_min = Networks.get_p_min(generator)

            for s in scenarios
                for ts in target_timepoints
                    # pmin B_on < P_tso_min < P_tso_max < pmax B_on
                    @constraint(tso_model, p_min * b_on_vars[gen_id, ts, s] <= p_imposition_min[gen_id, ts, s]);
                    @constraint(tso_model, p_imposition_min[gen_id, ts, s] <= p_imposition_max[gen_id, ts, s])
                    @constraint(tso_model, p_imposition_max[gen_id, ts, s] <= p_max * b_on_vars[gen_id, ts, s])
                end
            end
        end
    end
end

function add_tso_constraints!(bimodel_container::TSOBilevelModel,
                            target_timepoints, scenarios, network,
                            uncertainties_at_ech::UncertaintiesAtEch)
    tso_model_container::TSOBilevelTSOModelContainer = bimodel_container.upper
    market_model_container::TSOBilevelMarketModelContainer = bimodel_container.lower

    rso_constraints!(bimodel_container.model,
                    tso_model_container.flows,
                    tso_model_container.rso_constraint,
                    market_model_container.pilotable_model, #pilotable injections are decided by Market
                    tso_model_container.limitable_model, #limitable injections are decided by TSO
                    tso_model_container.lol_model,
                    target_timepoints, scenarios,
                    uncertainties_at_ech, network)
    add_loss_of_load_distribution_constraint!(tso_model_container,
                                        market_model_container.lol_model,
                                        target_timepoints, scenarios, network)
    add_enr_distribution_constraint!(tso_model_container,
                                    market_model_container.limitable_model,
                                    target_timepoints, scenarios, network)
    add_capping_distribution_constraint!(tso_model_container,
                                    market_model_container.limitable_model,
                                    target_timepoints, scenarios, network)
    add_injection_commitment_constraints!(tso_model_container,
                                    market_model_container.pilotable_model,
                                    target_timepoints, scenarios, network)

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
                            # generators_initial_state::SortedDict{String,GeneratorState},
                            scenarios::Vector{String},
                            uncertainties_at_ech::UncertaintiesAtEch,
                            firmness::Firmness,
                            # preceding_tso_schedule::Schedule,
                            # preceding_tso_actions::TSOActions,
                            configs::TSOBilevelConfigs
                            )
    add_limitables!(model_container,
                    network, target_timepoints, scenarios, uncertainties_at_ech)
    add_pilotables!(model_container,
                    network, target_timepoints, scenarios, firmness, configs.LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET)
    add_slacks!(model_container, network, target_timepoints, scenarios, uncertainties_at_ech)
end

function add_limitables!(model_container::TSOBilevelMarketModelContainer,
                        network::Networks.Network,
                        target_timepoints::Vector{Dates.DateTime},
                        scenarios::Vector{String},
                        uncertainties_at_ech::UncertaintiesAtEch
                        #firmness::Firmness,
                        #always_link_scenarios_limit=false
                        )
    model = model_container.model
    limitable_model = model_container.limitable_model
    for ts in target_timepoints
        for s in scenarios
            enr_max = compute_prod(uncertainties_at_ech, network, ts, s)
            name =  @sprintf("P_injected[%s,%s]", ts, s)
            limitable_model.p_injected[ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=enr_max)
            name =  @sprintf("P_capping[%s,%s]", ts, s)
            limitable_model.p_capping[ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=enr_max)
        end
    end
    return model_container
end

function add_pilotables!(model_container::TSOBilevelMarketModelContainer,
                        network::Networks.Network,
                        target_timepoints::Vector{Dates.DateTime},
                        scenarios::Vector{String},
                        #generators_initial_state::SortedDict{String,GeneratorState},
                        firmness::Firmness,
                        #reference_schedule::Schedule,
                        #tso_actions::TSOActions,
                        always_link_levels::Bool
                        )
    pilotable_generators = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    for pilotable_gen in pilotable_generators
        gen_id = Networks.get_id(pilotable_gen)
        add_pilotable!(model_container.pilotable_model, model_container.model,
                        pilotable_gen,
                        target_timepoints,
                        scenarios,
                        get_power_level_firmness(firmness, gen_id),
                        #no TSO actions, the upper problem considers them
                        always_link_levels
                        )
    end
    return model_container.pilotable_model
end
function add_pilotable!(pilotable_model::TSOBilevelMarketPilotableModel, model::Model,
                        generator::Networks.Generator,
                        target_timepoints,
                        scenarios,
                        power_level_firmness::SortedDict{Dates.DateTime, DecisionFirmness},
                        always_link_levels::Bool
                        )
    gen_id = Networks.get_id(generator)
    p_max = Networks.get_p_max(generator)
    for ts in target_timepoints
        for s in scenarios
            add_p_injected!(pilotable_model, model, gen_id, ts, s, p_max, false)
        end
    end

    add_scenarios_linking_constraints!(model, generator,
                                        pilotable_model.p_injected,
                                        target_timepoints, scenarios,
                                        power_level_firmness,
                                        always_link_levels
                                        )

    return pilotable_model, model
end

function add_slacks!(model_container::TSOBilevelMarketModelContainer,
                    network::Networks.Network,
                    target_timepoints::Vector{Dates.DateTime},
                    scenarios::Vector{String},
                    uncertainties_at_ech::UncertaintiesAtEch)
    model = model_container.model
    lol_model = model_container.lol_model
    for ts in target_timepoints
        for s in scenarios
            conso_max = compute_load(uncertainties_at_ech, network, ts, s)
            name =  @sprintf("P_loss_of_load[%s,%s]", ts, s)
            lol_model.p_global_loss_of_load[ts, s] = @variable(model, base_name=name,
                                                        lower_bound=0., upper_bound=conso_max)
        end
    end
    return lol_model
end

function add_firmness_duals(kkt_model_container::TSOBilevelKKTModelContainer,
                            target_timepoints, scenarios,
                            network,
                            firmness,
                            link_scenarios_pilotable_level_market)
    for gen_pilotable in Networks.get_generators_of_type(network, Networks.PILOTABLE)
        gen_id = Networks.get_id(gen_pilotable)
        for ts in target_timepoints
            decision_firmness_l = get_power_level_firmness(firmness, gen_id, ts)
            for s in scenarios
                name = @sprintf("c_firmness[%s,%s,%s]",gen_id,ts,s)

                if requires_linking(decision_firmness_l, link_scenarios_pilotable_level_market)
                    #create duals relative to firmness constraints of injected pilotable power in market
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

function add_link_capping_constraint!(market_model_container::TSOBilevelMarketModelContainer,
                                        kkt_model_container::TSOBilevelKKTModelContainer,
                                    tso_limitable_model::TSOBilevelTSOLimitableModel,
                                    target_timepoints, scenarios, network)
    market_model = market_model_container.model
    market_limitable_model = market_model_container.limitable_model
    for ts in target_timepoints
        for s in scenarios
            name = @sprintf("c_min_e[%s,%s]",ts,s)
            @constraint(market_model,
                        tso_limitable_model.p_global_capping[ts,s] <= market_limitable_model.p_capping[ts,s],
                        base_name=name)

            #create duals and indicators relative to TSO min capping constraint
            add_dual_and_indicator!(kkt_model_container.model,
                                    kkt_model_container.capping_duals, kkt_model_container.capping_indicators, (ts,s),
                                    name, true)
        end
    end
    return market_model_container
end

function add_link_loss_of_load_constraint!(market_model_container::TSOBilevelMarketModelContainer,
                                        kkt_model_container::TSOBilevelKKTModelContainer,
                                    tso_lol_model::TSOBilevelTSOLoLModel,
                                    target_timepoints, scenarios, network)
    market_model = market_model_container.model
    market_lol_model = market_model_container.lol_model
    for ts in target_timepoints
        for s in scenarios
            name = @sprintf("c_min_lol[%s,%s]",ts,s)
            @constraint(market_model,
                        tso_lol_model.p_global_loss_of_load[ts,s] <= market_lol_model.p_global_loss_of_load[ts,s],
                        base_name=name)

            #create duals and indicators relative to TSO min LoL constraint
            add_dual_and_indicator!(kkt_model_container.model,
                                    kkt_model_container.loss_of_load_duals, kkt_model_container.loss_of_load_indicators, (ts,s),
                                    name, true)
        end
    end
    return market_model_container
end
function add_pilotables_constraints!(market_model_container::TSOBilevelMarketModelContainer,
                                    kkt_model_container::TSOBilevelKKTModelContainer,
                                tso_pilotable_model::TSOBilevelTSOPilotableModel,
                                target_timepoints, scenarios, network)
    for ts in target_timepoints
        for s in scenarios
            for gen in Networks.get_generators_of_type(network, Networks.PILOTABLE)
                gen_id = Networks.get_id(gen)
                add_pilotable_pmin_constraints!(market_model_container, kkt_model_container,
                                            tso_pilotable_model, gen_id, ts, s)
                add_pilotable_pmax_constraints!(market_model_container, kkt_model_container,
                                            tso_pilotable_model, gen_id, ts, s)
            end
        end
    end

    return market_model_container, kkt_model_container
end
function add_pilotable_pmin_constraints!(market_model_container::TSOBilevelMarketModelContainer,
                                        kkt_model_container::TSOBilevelKKTModelContainer,
                                        tso_pilotable_model::TSOBilevelTSOPilotableModel,
                                        gen_id, ts, s)
    injection_var = market_model_container.pilotable_model.p_injected[gen_id,ts,s]
    name = @sprintf("c_tso_pmin[%s,%s,%s]",gen_id,ts,s)
    @constraint(market_model_container.model,
                tso_pilotable_model.p_imposition_min[gen_id,ts,s] <= injection_var,
                base_name=name)
    #create duals and indicators relative to tso pmin constraint
    add_dual_and_indicator!(kkt_model_container.model,
                            kkt_model_container.pmin_duals, kkt_model_container.pmin_indicators, (gen_id,ts,s),
                            name, true)

    return kkt_model_container
end
function add_pilotable_pmax_constraints!(market_model_container::TSOBilevelMarketModelContainer,
                                        kkt_model_container::TSOBilevelKKTModelContainer,
                                        tso_pilotable_model::TSOBilevelTSOPilotableModel,
                                        gen_id, ts, s)
    injection_var = market_model_container.pilotable_model.p_injected[gen_id,ts,s]
    name = @sprintf("c_tso_pmax[%s,%s,%s]",gen_id,ts,s)
    @constraint(market_model_container.model,
                injection_var <= tso_pilotable_model.p_imposition_max[gen_id,ts,s],
                base_name=name)
    #create duals and indicators relative to tso pmax constraint
    add_dual_and_indicator!(kkt_model_container.model,
                            kkt_model_container.pmax_duals, kkt_model_container.pmax_indicators, (gen_id,ts,s),
                            name, true)

    return kkt_model_container
end

function add_market_constraints!(bimodel_container::TSOBilevelModel,
                            target_timepoints, scenarios, network,
                            uncertainties_at_ech::UncertaintiesAtEch)
    tso_model_container::TSOBilevelTSOModelContainer = bimodel_container.upper
    market_model_container::TSOBilevelMarketModelContainer = bimodel_container.lower
    kkt_model_container::TSOBilevelKKTModelContainer = bimodel_container.kkt_model

    eod_constraints!(market_model_container.model, market_model_container.eod_constraint,
                    market_model_container.pilotable_model,
                    tso_model_container.limitable_model,
                    market_model_container.lol_model,
                    target_timepoints, scenarios,
                    uncertainties_at_ech, network
                    )
    add_eod_constraints_duals!(kkt_model_container, target_timepoints, scenarios)

    add_link_capping_constraint!(market_model_container, kkt_model_container,
                                tso_model_container.limitable_model,
                                target_timepoints, scenarios, network)
    add_link_loss_of_load_constraint!(market_model_container, kkt_model_container,
                                tso_model_container.lol_model,
                                target_timepoints, scenarios, network)
    add_pilotables_constraints!(market_model_container, kkt_model_container,
                            tso_model_container.pilotable_model,
                            target_timepoints, scenarios, network)

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
    objective_model.limitable_cost += coeffxsum(model_container.limitable_model.p_capping, capping_cost)

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

function add_dual_and_indicator!(kkt_model::Model, duals_dict, indicators_dict,
                                key, name,
                                is_positive::Bool)

    add_dual!(kkt_model, duals_dict, key, name, is_positive)

    indicators_dict[key] = @variable(kkt_model, binary=true, base_name="indicator_"*name)

    return duals_dict[key] , indicators_dict[key]
end

function add_dual!(kkt_model::Model, duals_dict,
                    key, name,
                    is_positive::Bool)
    if is_positive
        duals_dict[key] = @variable(kkt_model, lower_bound=0., base_name="dual_"*name)
    else
        duals_dict[key] = @variable(kkt_model, base_name="dual_"*name)
    end

    return duals_dict[key]
end

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
            #                                                 market_model.limitable_model.p_capping[ts,s]) )
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

function add_emin_complementarity_constraints!(model_container::TSOBilevelModel,
                                            big_m, target_timepoints, scenarios)
    tso_limitable_model = model_container.upper.limitable_model
    market_limitable_model = model_container.lower.limitable_model
    kkt_model = model_container.kkt_model

    for ts in target_timepoints
        for s in scenarios
            kkt_var = kkt_model.capping_duals[ts,s]
            cstr_expr = market_limitable_model.p_capping[ts,s] - tso_limitable_model.p_global_capping[ts,s]
            b_indicator = kkt_model.capping_indicators[ts,s]
            ub_cstr = compute_ub(cstr_expr, big_m)
            formulate_complementarity_constraints!(kkt_model.model, kkt_var, cstr_expr, b_indicator, big_m, ub_cstr)
        end
    end
    return model_container
end

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

    bimodel_container_l = TSOBilevelModel()
    @assert(configs.big_m >= configs.MARKET_LOL_PENALTY)
    @assert(configs.big_m >= configs.MARKET_CAPPING_COST)
    all( configs.big_m >= Networks.get_prop_cost(gen)
        for gen in Networks.get_generators_of_type(network, Networks.PILOTABLE) )

    if is_market(configs.REF_SCHEDULE_TYPE_IN_TSO)
        reference_schedule = preceding_market_schedule
    elseif is_tso(configs.REF_SCHEDULE_TYPE_IN_TSO)
        reference_schedule = preceding_tso_schedule
    else
        throw( error("Invalid REF_SCHEDULE_TYPE_IN_TSO config.") )
    end
    create_tso_vars!(bimodel_container_l.upper,
                    network, target_timepoints, generators_initial_state, scenarios,
                    uncertainties_at_ech, firmness, reference_schedule,
                    configs)
    create_market_vars!(bimodel_container_l.lower,
                        network, target_timepoints, scenarios, uncertainties_at_ech, firmness, configs)

    #this is the expression no objective is added to the jump model
    create_market_objectives!(bimodel_container_l.lower, network,
                            configs.MARKET_CAPPING_COST, configs.MARKET_LOL_PENALTY)

    #constraints may use upper and lower vars at the same time
    add_tso_constraints!(bimodel_container_l, target_timepoints, scenarios, network, uncertainties_at_ech)
    #kkt primal feasibility + variables creation
    add_firmness_duals(bimodel_container_l.kkt_model,
                        target_timepoints, scenarios, network, firmness,
                        configs.LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET)
    add_market_constraints!(bimodel_container_l, target_timepoints, scenarios, network, uncertainties_at_ech)
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

    solve!(bimodel_container_l, configs.problem_name, configs.out_path)
    @info("Lower Objective Value : $(value(bimodel_container_l.lower.objective_model.full_obj))")

    return bimodel_container_l
end
