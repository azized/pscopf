using .Networks

using JuMP
using Dates
using DataStructures
using Printf
using Parameters


"""
REF_SCHEDULE_TYPE : Indicates wether to consider the preceding market or TSO schedule as a reference.
                    The reference schedule is used to get decided commitment and production levels if
                      tso actions are missing.
"""
@with_kw mutable struct TSOConfigs <: AbstractRunnableConfigs
    CONSIDER_N_1_CSTRS::Bool = get_config("CONSIDER_N_1")
    loss_of_load_penalty = get_config("tso_loss_of_load_penalty_value")
    limitation_penalty = get_config("tso_limit_penalty_value")
    out_path = nothing
    problem_name = "TSO"
    REF_SCHEDULE_TYPE::Union{Market,TSO} = TSO();
end


@with_kw struct TSOLimitableModel <: AbstractLimitableModel
    #gen,ts,s
    delta_p = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_injected = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_limit = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_is_limited = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_limit_x_is_limited = SortedDict{Tuple{String,DateTime,String},VariableRef}();
end

@with_kw struct TSOPilotableModel <: AbstractPilotableModel
    #gen,ts,s
    delta_p = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    p_injected = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_start = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_on = SortedDict{Tuple{String,DateTime,String},VariableRef}();
end

@with_kw struct TSOLoLModel <: AbstractLoLModel
    #bus,ts,s
    p_loss_of_load = SortedDict{Tuple{String,DateTime,String},VariableRef}();
end

@with_kw mutable struct TSOObjectiveModel <: AbstractObjectiveModel
    deltas = AffExpr(0)

    prop_cost = AffExpr(0)
    start_cost = AffExpr(0)

    penalty = AffExpr(0)

    full_obj_1 = AffExpr(0)
    full_obj_2 = AffExpr(0)
end


@with_kw mutable struct TSOModel <: AbstractModelContainer
    model::Model = Model()
    limitable_model::TSOLimitableModel = TSOLimitableModel()
    pilotable_model::TSOPilotableModel = TSOPilotableModel()
    lol_model::TSOLoLModel = TSOLoLModel()
    objective_model::TSOObjectiveModel = TSOObjectiveModel()
    #ts,s
    eod_constraint::SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef} =
        SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef}()
    #branch,ts,s,ptdf_case
    rso_combinations::Vector{Tuple{String,DateTime,String,String}} =
        Vector{Tuple{String,DateTime,String,String}}()
    flows::SortedDict{Tuple{String,DateTime,String,String},AffExpr} =
        SortedDict{Tuple{String,DateTime,String,String},AffExpr}()
    rso_constraints::SortedDict{Tuple{String,DateTime,String,String},Tuple{ConstraintRef,ConstraintRef}} =
        SortedDict{Tuple{String,DateTime,String,String},Tuple{ConstraintRef,ConstraintRef}}()
    #deltas constraint
    deltas_bounding_constraint::Union{ConstraintRef, Missing} = missing
end

function has_positive_slack(model_container::TSOModel)::Bool
    return has_positive_value(model_container.lol_model.p_loss_of_load)
end


function sum_capping(limitable_model::TSOLimitableModel, ts,s, network::Networks.Network)
    error("TODO : requires uncertainties cause capping=uncertainties-injection")
end

#TODO define a struct LocalisedLolModel to use it for TSOLoLModel and TSOBilevelTSOLoLModel
function sum_lol(lol_model::TSOLoLModel, ts, s, network::Networks.Network)
    sum_l = 0.
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        sum_l += lol_model.p_loss_of_load[bus_id,ts,s]
    end
    return sum_l
end


function add_tso_limitable_prop_cost!(obj_component::AffExpr,
                                uncertainties_at_ech::UncertaintiesAtEch,
                                p_injected::AbstractDict{T,V}, network)  where T <: Tuple where V <: VariableRef
    #NOTE: need to make sure uncertainty > injection

    for ((gen_id,ts,s), p_injected_var) in p_injected
        generator = Networks.get_generator(network, gen_id)
        gen_prop_cost = Networks.get_prop_cost(generator)
        uncertainty = get_uncertainties(uncertainties_at_ech, gen_id, ts, s)
        add_to_expression!(obj_component,
                            (uncertainty - p_injected_var) * gen_prop_cost)
    end

    return obj_component
end

function create_objectives!(model_container::TSOModel,
                            network, uncertainties_at_ech, gratis_starts, loss_of_load_cost, limitation_penalty)

    # cost for cutting load/consumption
    add_coeffxsum_cost!(model_container.objective_model.penalty,
                        model_container.lol_model.p_loss_of_load, loss_of_load_cost)

    # avoid limiting when not necessary
    for (_, var_is_limited) in model_container.limitable_model.b_is_limited
        model_container.objective_model.penalty += limitation_penalty * var_is_limited
    end

    ## Objective 1 :

    # cost for deviating from market schedule
    for (_, var_delta) in model_container.pilotable_model.delta_p
        model_container.objective_model.deltas += var_delta
    end
    for (_, var_delta) in model_container.limitable_model.delta_p
        model_container.objective_model.deltas += var_delta
    end

    ## Objective 2 :

    # cost for starting pilotables
    add_pilotable_start_cost!(model_container.objective_model.start_cost,
                            model_container.pilotable_model.b_start, network, gratis_starts)

    # cost for limitables : cost of capped limitable power
    add_tso_limitable_prop_cost!(model_container.objective_model.prop_cost,
                                uncertainties_at_ech,
                                model_container.limitable_model.p_injected, network)

    # cost for using pilotables
    add_prop_cost!(model_container.objective_model.prop_cost,
                            model_container.pilotable_model.p_injected, network)

    # Objective 1 :
    model_container.objective_model.full_obj_1 = ( model_container.objective_model.deltas +
                                                model_container.objective_model.penalty )
    # Objective 2 :
    model_container.objective_model.full_obj_2 = ( model_container.objective_model.start_cost +
                                                model_container.objective_model.prop_cost +
                                                model_container.objective_model.penalty )

    return model_container
end

function create_tso_model(network::Networks.Network,
                        target_timepoints::Vector{Dates.DateTime},
                        generators_initial_state::SortedDict{String,GeneratorState},
                        scenarios::Vector{String},
                        uncertainties_at_ech::UncertaintiesAtEch,
                        firmness::Firmness,
                        preceding_market_schedule::Schedule,
                        preceding_tso_schedule::Schedule,
                        gratis_starts::Set{Tuple{String,Dates.DateTime}},
                        configs::TSOConfigs)
    # TODO : check coherence between : preceding_reference_schedule and TSOActions.impositions cause we do not consider TSOActions
    if is_market(configs.REF_SCHEDULE_TYPE)
        reference_schedule = preceding_market_schedule
    elseif is_tso(configs.REF_SCHEDULE_TYPE)
        reference_schedule = preceding_tso_schedule
    else
        throw( error("Invalid REF_SCHEDULE_TYPE config.") )
    end

    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    limitables_list_l = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    buses_list_l = Networks.get_buses(network)

    model_container_l = TSOModel()

    # Variables
    add_pilotables_vars!(model_container_l,
                            pilotables_list_l, target_timepoints, scenarios,
                            preceding_market_schedule,
                            injection_vars=true, commitment_vars=true, delta_vars=true)
    add_limitables_vars!(model_container_l, target_timepoints, scenarios,
                            limitables_list_l, preceding_market_schedule,
                            injection_vars=true, limit_vars=true, delta_vars=true
                            )
    add_lol_vars!(model_container_l, target_timepoints, scenarios,
                    buses_list_l,
                    local_lol_vars=true)


    # Constraints

    # Pilotables
    pilotable_power_constraints!(model_container_l.model,
                                model_container_l.pilotable_model, pilotables_list_l, target_timepoints, scenarios,
                                firmness, reference_schedule,
                                always_link_scenarios=false)
    unit_commitment_constraints!(model_container_l.model,
                                model_container_l.pilotable_model, pilotables_list_l,  target_timepoints, scenarios,
                                firmness, reference_schedule, generators_initial_state,
                                always_link_scenarios=false)

    # Limitables
    limitable_power_constraints!(model_container_l.model,
                                model_container_l.limitable_model, limitables_list_l, target_timepoints, scenarios,
                                firmness, uncertainties_at_ech,
                                always_link_scenarios=false)

    # LoL
    local_lol_constraints!(model_container_l.model,
                            model_container_l.lol_model, buses_list_l, target_timepoints, scenarios,
                            uncertainties_at_ech)

    # EOD
    eod_constraints!(model_container_l.model, model_container_l.eod_constraint,
                    model_container_l.pilotable_model,
                    model_container_l.limitable_model,
                    model_container_l.lol_model,
                    target_timepoints, scenarios,
                    uncertainties_at_ech, network
                    )

    # RSO
    #just adds the combinations to consider not the constraints
    add_rso_combinations!(get_rso_combinations(model_container_l),
                        configs.CONSIDER_N_1_CSTRS,
                        network, target_timepoints, scenarios)

    create_objectives!(model_container_l,
                        network, uncertainties_at_ech,
                        gratis_starts,
                        configs.loss_of_load_penalty, configs.limitation_penalty)

    return model_container_l
end

function tso_out_fo(network::Networks.Network,
                    target_timepoints::Vector{Dates.DateTime},
                    generators_initial_state::SortedDict{String,GeneratorState},
                    scenarios::Vector{String},
                    uncertainties_at_ech::UncertaintiesAtEch,
                    firmness::Firmness,
                    preceding_market_schedule::Schedule,
                    preceding_tso_schedule::Schedule,
                    gratis_starts::Set{Tuple{String,Dates.DateTime}},
                    configs::TSOConfigs
                    )
    @assert all(configs.loss_of_load_penalty > Networks.get_prop_cost(gen)
                for gen in Networks.get_generators(network))
    @assert all(configs.loss_of_load_penalty > Networks.get_prop_cost(gen) + Networks.get_start_cost(gen)/Networks.get_p_min(gen)
                for gen in Networks.get_generators(network)
                if Networks.needs_commitment(gen))

    @timeit TIMER_TRACKS "tso_modeling" model_container_l = create_tso_model(network, target_timepoints, generators_initial_state,
                                                                            scenarios, uncertainties_at_ech, firmness,
                                                                            preceding_market_schedule, preceding_tso_schedule,
                                                                            gratis_starts,
                                                                            configs)

    @timeit TIMER_TRACKS "tso_solve" tso_solve!(model_container_l,
                                            solve_2steps_deltas!, configs,
                                            uncertainties_at_ech, network,
                                            get_config("ADD_RSO_CSTR_DYNAMICALLY"))
    log_flows(model_container_l, network, configs.out_path, configs.problem_name)

    return model_container_l
end
