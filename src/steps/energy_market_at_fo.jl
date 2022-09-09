using .Networks

using Dates
using JuMP
using Statistics
using Printf

@with_kw mutable struct EnergyMarketAtFO <: AbstractMarket
    configs = EnergyMarketConfigs()
end

function aggregate_scenarios(context::AbstractContext, ech::Dates.DateTime)
    agg_scenario_name = aggregate_scenario_name(context, ech)

    agg_uncertainties_at_ech = UncertaintiesAtEch()
    uncertainties_at_ech = get_uncertainties(context, ech)
    for (injection_name, _) in uncertainties_at_ech
        for (ts, by_scenario_injections) in get_uncertainties(uncertainties_at_ech, injection_name)
            values = []
            for (_, val) in by_scenario_injections
                push!(values, val)
            end
            value_l = mean(values)
            PSCOPF.add_uncertainty!(agg_uncertainties_at_ech, injection_name, ts, agg_scenario_name, value_l)
        end
    end

    return agg_scenario_name, agg_uncertainties_at_ech
end

function run(runnable::EnergyMarketAtFO,
            ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime},
            context::AbstractContext)
    fo_start_time = TS[1] - get_fo_length(get_management_mode(context))
    if fo_start_time != ech
        msg = @sprintf("invalid step at ech=%s : EnergyMarketAtFO needs to be launched at FO start (ie %s)", ech, fo_start_time)
        @warn(msg)
    end

    problem_name_l = @sprintf("energy_market_at_FO_%s", ech)

    tso_actions = filter_tso_actions(get_tso_actions(context))

    gratis_starts = Set{Tuple{String,Dates.DateTime}}()
    if runnable.configs.CONSIDER_GRATIS_STARTS
        gratis_starts = init_gratis_start(context, runnable.configs.REF_SCHEDULE_TYPE)
    end

    agg_scenario_name, agg_uncertainties = aggregate_scenarios(context, ech)

    @assert(length(get_scenarios(agg_uncertainties)) == 1)
    @assert check_uncertainties(Uncertainties(ech=>agg_uncertainties), get_network(context))
    @assert check_uncertainties_contain_ts(Uncertainties(ech=>agg_uncertainties), get_target_timepoints(context))

    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    return energy_market(get_network(context),
                        TS,
                        get_generators_initial_state(context),
                        [agg_scenario_name],
                        agg_uncertainties,
                        firmness,
                        get_market_schedule(context), #this uses original scenario names
                        get_tso_schedule(context),
                        tso_actions,
                        gratis_starts,
                        runnable.configs
                        )
end

function update_market_schedule!(context::AbstractContext, ech,
                                result::EnergyMarketModel,
                                firmness,
                                runnable::EnergyMarketAtFO)
    market_schedule = get_market_schedule(context)

    market_schedule.decider_type = DeciderType(runnable)
    market_schedule.decision_time = ech

    for ((gen_id, ts, s), p_injected_var) in result.pilotable_model.p_injected
        set_prod_definitive_value!(market_schedule, gen_id, ts, value(p_injected_var))
        if get_power_level_firmness(firmness, gen_id, ts) == DECIDED
            @assert( ismissing(get_prod_value(market_schedule, gen_id, ts))
                    || value(p_injected_var) ≈ get_prod_value(market_schedule, gen_id, ts) )
            set_prod_definitive_value!(market_schedule, gen_id, ts, value(p_injected_var))
        end
    end

    for ((gen_id, ts, _), b_on_var) in result.pilotable_model.b_on
        gen_state_value = parse(GeneratorState, value(b_on_var))
        set_commitment_definitive_value!(market_schedule, gen_id, ts, gen_state_value)
    end

    # TODO : may need to adapt cause result handles only one scenario
    # Capping
    update_schedule_capping!(market_schedule, context, ech, result.limitable_model, runnable.configs.CONSIDER_TSOACTIONS_LIMITATIONS)

    # Limitables levels : needs to be after capping update
    #FIXME : avoid re-aggregating scenarios
    agg_scenario_name, agg_uncertainties = aggregate_scenarios(context, ech)
    for gen in Networks.get_generators_of_type(get_network(context), Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        for ts in get_target_timepoints(context)
            for s in get_scenarios(context)
                capped_l = safeget_capping(market_schedule, gen_id, ts, s)
                injected_l = get_uncertainties(agg_uncertainties,gen_id,ts,agg_scenario_name) - capped_l
                set_prod_value!(market_schedule, gen_id, ts, s, injected_l)
            end
        end
    end

    # loss_of_load (load-shedding)
    update_schedule_loss_of_load!(market_schedule, context, ech, result.lol_model)

    return market_schedule
end
