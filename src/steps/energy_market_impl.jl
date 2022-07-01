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
@with_kw mutable struct EnergyMarketConfigs <: AbstractRunnableConfigs
    loss_of_load_penalty = get_config("market_loss_of_load_penalty_value")
    out_path = nothing
    problem_name = "EnergyMarket"
    REF_SCHEDULE_TYPE::Union{Market,TSO} = TSO();
    CONSIDER_TSOACTIONS_LIMITATIONS::Bool = false
    CONSIDER_TSOACTIONS_IMPOSITIONS::Bool = false
    CONSIDER_GRATIS_STARTS::Bool = true
end

@with_kw struct EnergyMarketLimitableModel <: AbstractLimitableModel
    #ts,s
    p_global_capping = SortedDict{Tuple{DateTime,String},VariableRef}();
    #firmness_constraints
end

@with_kw struct EnergyMarketPilotableModel <: AbstractPilotableModel
    #gen,ts,s
    p_injected = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_start = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #gen,ts,s
    b_on = SortedDict{Tuple{String,DateTime,String},VariableRef}();
    #commitment_constraints = Dict{Tuple{String,DateTime,String},ConstraintRef}();
    #firmness_constraints
end

@with_kw struct EnergyMarketLoLModel <: AbstractLoLModel
    #ts,s
    p_global_loss_of_load = SortedDict{Tuple{DateTime,String},VariableRef}();
end

@with_kw mutable struct EnergyMarketObjectiveModel <: AbstractObjectiveModel
    prop_cost = AffExpr(0)
    start_cost = AffExpr(0)
    penalty = AffExpr(0)

    full_obj = AffExpr(0)
end


@with_kw mutable struct EnergyMarketModel <: AbstractModelContainer
    model::Model = Model()
    limitable_model::EnergyMarketLimitableModel = EnergyMarketLimitableModel()
    pilotable_model::EnergyMarketPilotableModel = EnergyMarketPilotableModel()
    lol_model::EnergyMarketLoLModel = EnergyMarketLoLModel()
    objective_model::EnergyMarketObjectiveModel = EnergyMarketObjectiveModel()
    eod_constraint::SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef} =
        SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef}()
end

SCENARIOS_DELIMITER = "_+_"

function aggregate_scenario_name(scenarios::Vector{String})
    return join(scenarios, SCENARIOS_DELIMITER)
end
function aggregate_scenario_name(context::AbstractContext, ech::Dates.DateTime)
    scenarios = get_scenarios(context, ech)
    return aggregate_scenario_name(scenarios)
end

function has_positive_slack(model_container::EnergyMarketModel)::Bool
    return has_positive_value(model_container.lol_model.p_global_loss_of_load)
end

"""
    energy_market
# Arguments
    - `network::Networks.Network`
    - `target_timepoints::Vector{Dates.DateTime}`
    - `generators_initial_state::SortedDict{String,GeneratorState}` :
        The ON/OFF state of each generator just before the first target timepoint.
    - `scenarios::Vector{String}` : considered scenarios
    - `uncertainties_at_ech::UncertaintiesAtEch` :
        The considered limitables injections and bus consumption realisations per scenario
         for the current timepoint of execution (i.e. `ech`).
    - `firmness::Firmness` : The required level of firmness for commitment and power level decisions
    - `preceding_market_schedule::Schedule` : The preceding market schedule may be used to set already decided values.
    - `preceding_tso_schedule::Schedule` : The preceding tso schedule may be used to set already decided values.
    - `tso_actions::TSOActions` : The preceding tso actions maybe used for additional constraints.
        The behaviour of the model wit hrespect to the TSO actions is defined by the parameters
         `CONSIDER_TSOACTIONS_LIMITATIONS` and `CONSIDER_TSOACTIONS_IMPOSITIONS`
         of `configs::EnergyMarketConfigs`
    - `gratis_starts::Set{Tuple{String,Dates.DateTime}}` :
        Tuples of (gen_id, ts) giving already paid for starting decisions. So if the market starts
         unit gen_id at timestep ts, it will not pay the starting cost.
        These are ignored if EnergyMarketConfigs::CONSIDER_GRATIS_STARTS is false
    - `configs::EnergyMarketConfigs` :
        settings modifying the behaviour of the model.
"""
function energy_market(network::Networks.Network,
                    target_timepoints::Vector{Dates.DateTime},
                    generators_initial_state::SortedDict{String,GeneratorState},
                    scenarios::Vector{String},
                    uncertainties_at_ech::UncertaintiesAtEch,
                    firmness::Firmness,
                    preceding_market_schedule::Schedule,
                    preceding_tso_schedule::Schedule,
                    tso_actions::TSOActions,
                    gratis_starts::Set{Tuple{String,Dates.DateTime}},
                    configs::EnergyMarketConfigs
                    )
    #if gratis_starts is not empty, we should CONSIDER_GRATIS_STARTS
    @assert( configs.CONSIDER_GRATIS_STARTS | isempty(gratis_starts) )
    #if not there a risk of cutting conso instead of using a generator
    @assert all(configs.loss_of_load_penalty > Networks.get_prop_cost(gen)
                for gen in Networks.get_generators(network))
    @assert all(configs.loss_of_load_penalty > Networks.get_prop_cost(gen) + Networks.get_start_cost(gen)/Networks.get_p_min(gen)
                    for gen in Networks.get_generators(network)
                    if Networks.needs_commitment(gen))

    @timeit TIMER_TRACKS "market_modeling" model_container_l = create_market_model(network,
                                                                                target_timepoints, generators_initial_state,
                                                                                scenarios, uncertainties_at_ech, firmness,
                                                                                preceding_market_schedule, preceding_tso_schedule,
                                                                                tso_actions, gratis_starts, configs)

    @timeit TIMER_TRACKS "market_solve" solve!(model_container_l, configs.problem_name, configs.out_path)

    return model_container_l
end

function create_market_model(network::Networks.Network,
                            target_timepoints::Vector{Dates.DateTime},
                            generators_initial_state::SortedDict{String,GeneratorState},
                            scenarios::Vector{String},
                            uncertainties_at_ech::UncertaintiesAtEch,
                            firmness::Firmness,
                            preceding_market_schedule::Schedule,
                            preceding_tso_schedule::Schedule,
                            tso_actions::TSOActions,
                            gratis_starts::Set{Tuple{String,Dates.DateTime}},
                            configs::EnergyMarketConfigs)
    if is_market(configs.REF_SCHEDULE_TYPE)
        reference_schedule = preceding_market_schedule
    elseif is_tso(configs.REF_SCHEDULE_TYPE)
        reference_schedule = preceding_tso_schedule
    else
        throw( error("Invalid REF_SCHEDULE_TYPE config.") )
    end

    pilotables_list_l = Networks.get_generators_of_type(network, Networks.PILOTABLE)
    limitables_list_l = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    buses_list = Networks.get_buses(network)

    model_container_l = EnergyMarketModel()

    # Variables
    add_pilotables_vars!(model_container_l,
                        pilotables_list_l, target_timepoints, scenarios,
                        injection_vars=true, commitment_vars=true)
    add_limitables_vars!(model_container_l, target_timepoints, scenarios,
                        global_capping_vars=true
                        )
    add_lol_vars!(model_container_l, target_timepoints, scenarios,
                global_lol_vars=true)

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
    if configs.CONSIDER_TSOACTIONS_IMPOSITIONS
        respect_impositions_constraints!(model_container_l.model,
                                        model_container_l.pilotable_model, pilotables_list_l,  target_timepoints, scenarios,
                                        tso_actions)
    end

    # Limitables
    if configs.CONSIDER_TSOACTIONS_LIMITATIONS
        global_capping_constraints!(model_container_l.model,
                                    model_container_l.limitable_model,limitables_list_l, target_timepoints, scenarios,
                                    uncertainties_at_ech, tso_actions=tso_actions)
    else
        global_capping_constraints!(model_container_l.model,
                                    model_container_l.limitable_model, limitables_list_l, target_timepoints, scenarios,
                                    uncertainties_at_ech)
    end

    # LoL
    global_lol_constraints!(model_container_l.model,
                            model_container_l.lol_model, buses_list, target_timepoints, scenarios,
                            uncertainties_at_ech)

    # EOD
    eod_constraints!(model_container_l.model, model_container_l.eod_constraint,
                        model_container_l.pilotable_model,
                        model_container_l.limitable_model,
                        model_container_l.lol_model,
                        target_timepoints, scenarios,
                        uncertainties_at_ech, network)


    add_objective!(model_container_l, network, gratis_starts, configs.loss_of_load_penalty)

    return model_container_l
end

function add_objective!(model_container::EnergyMarketModel, network, gratis_starts, loss_of_load_cost)
    # cost for starting pilotables
    add_pilotable_start_cost!(model_container.objective_model.start_cost,
                            model_container.pilotable_model.b_start, network, gratis_starts)

    # No limitable cost

    # cost for using pilotables
    add_prop_cost!(model_container.objective_model.prop_cost,
                    model_container.pilotable_model.p_injected, network)

    # cost for cutting load/consumption
    add_coeffxsum_cost!(model_container.objective_model.penalty,
                        model_container.lol_model.p_global_loss_of_load, loss_of_load_cost)

    model_container.objective_model.full_obj = ( model_container.objective_model.start_cost +
                                                model_container.objective_model.prop_cost +
                                                model_container.objective_model.penalty )
    @objective(model_container.model, Min, model_container.objective_model.full_obj)
    return model_container
end


###########################
# Context-update
###########################

function update_schedule_capping!(market_schedule, context, ech, limitable_model::EnergyMarketLimitableModel, consider_limitations::Bool)
    reset_capping!(market_schedule)
    uncertainties = get_uncertainties(context, ech)

    for ((ts, scenario), p_capping_var) in limitable_model.p_global_capping
        capped_power = value(p_capping_var)
        for s in split_str(scenario, SCENARIOS_DELIMITER, keepempty=false)
        #for EnergyMarket, s==scenario
        #for EnergyMarketAtFO, this allows handling aggregate scenarios "S1_+_S2"
            limitables_ids = Networks.get_id.(Networks.get_generators_of_type(get_network(context), Networks.LIMITABLE))

            if capped_power > 1e-09
            #distribute the capped power on limitables
                limitations = consider_limitations ? get_limitations(get_tso_actions(context)) : Limitations()
                capped_by_limitations = compute_capped(uncertainties, limitations, get_network(context), ts, s)
                capped_by_eod = capped_power - capped_by_limitations

                @printf("capped %f power in scenario %s at ts %s\n", capped_power, s, ts)
                @printf("capped %f power for eod reasons in scenario %s at ts %s\n", capped_by_eod, s, ts)
                distribution_key = Dict{String,Float64}(gen_id_l => get_capacity(gen_id_l, ts, s,
                                                                                limitations,
                                                                                uncertainties)
                                                        for gen_id_l in limitables_ids)
                distribution_key = normalize_values(distribution_key)
                @debug(@sprintf("distribution_key : %s", distribution_key))

                for (gen_id, coeff) in distribution_key
                    capped_value = ( get_capped_by_limitations(gen_id, ts, s, limitations, uncertainties)
                                    + coeff * capped_by_eod )
                    market_schedule.capping[gen_id, ts, s] = capped_value
                end
            else
                for gen_id in limitables_ids
                    market_schedule.capping[gen_id, ts, s] = 0.
                end
            end
        end
    end
end

function update_schedule_loss_of_load!(market_schedule, context, ech, lol_model::EnergyMarketLoLModel)
    reset_loss_of_load_by_bus!(market_schedule)
    uncertainties = get_uncertainties(context, ech)

    for ((ts, scenario), p_loss_of_load_var) in lol_model.p_global_loss_of_load
        total_loss_of_load = value(p_loss_of_load_var)
        for s in split_str(scenario, SCENARIOS_DELIMITER, keepempty=false)
        #for EnergyMarket, s==scenario
        #for EnergyMarketAtFO, this allows handling aggregate scenarios "S1_+_S2"
            bus_ids = Networks.get_id.(Networks.get_buses(get_network(context)))

            if total_loss_of_load > 1e-09
            #distribute the cut conso on buses
                @printf("cut conso %f in scenario %s at ts %s\n", total_loss_of_load, s, ts)
                distribution_key = Dict{String,Float64}(bus_id_l => get_uncertainties(uncertainties, bus_id_l, ts, s)
                                                        for bus_id_l in bus_ids)
                distribution_key = normalize_values(distribution_key)
                @debug(@sprintf("distribution_key : %s", distribution_key))

                for (bus_id, coeff) in distribution_key
                    cut_load_on_bus = coeff * total_loss_of_load
                    market_schedule.loss_of_load_by_bus[bus_id, ts, s] = cut_load_on_bus
                end
            else
                for bus_id in bus_ids
                    market_schedule.loss_of_load_by_bus[bus_id, ts, s] = 0.
                end
            end
        end
    end
end

function get_capacity(gen_id, ts, s, limitations, uncertainties_at_ech)
    p_lim = get_limitation(limitations, gen_id, ts, s)

    if ismissing(p_lim)
        return get_uncertainties(uncertainties_at_ech, gen_id, ts, s)
    else
        return min(p_lim, get_uncertainties(uncertainties_at_ech, gen_id, ts, s))
    end
end

function get_capped_by_limitations(gen_id, ts, s, limitations, uncertainties_at_ech)
    p_lim = get_limitation(limitations, gen_id, ts, s)

    gen_capped = 0.
    if !ismissing(p_lim) && get_uncertainties(uncertainties_at_ech, gen_id, ts, s) > p_lim
        gen_capped = ( get_uncertainties(uncertainties_at_ech, gen_id, ts, s) - p_lim )
        @debug(@sprintf("%s capped %f : limit is %f out of %f",
                        gen_id, gen_capped, p_lim, get_uncertainties(uncertainties_at_ech, gen_id, ts, s)))
    end
    return gen_capped
end


function compute_capped(uncertainties_at_ech::UncertaintiesAtEch,
                        limitations::SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Float64}},
                        limitable_generators, ts, s)
    if isempty(limitations)
        return 0.
    end

    capped = 0.
    for limitable_gen in limitable_generators
        gen_id = Networks.get_id(limitable_gen)
        capped += get_capped_by_limitations(gen_id, ts, s, limitations, uncertainties_at_ech)
    end
    return capped
end
function compute_capped(uncertainties_at_ech::UncertaintiesAtEch,
                        limitations::SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Float64}},
                        network::Networks.Network, ts, s)
    limitable_generators = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    return compute_capped(uncertainties_at_ech, limitations, limitable_generators, ts, s)
end
