using .Networks

using Dates
using JuMP
using Printf

@with_kw mutable struct EnergyMarket <: AbstractMarket
    configs = EnergyMarketConfigs()
end

function run(runnable::Union{EnergyMarket,BalanceMarket},
            ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime},
            context::AbstractContext)
    problem_name_l = @sprintf("%s_%s", typeof(runnable), ech)

    gratis_starts = Set{Tuple{String,Dates.DateTime}}()
    if runnable.configs.CONSIDER_GRATIS_STARTS
        gratis_starts = init_gratis_start(context, runnable.configs.REF_SCHEDULE_TYPE)
    end
    @debug("gratis_starts : ", gratis_starts)

    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    return energy_market(get_network(context),
                        TS,
                        get_generators_initial_state(context),
                        get_scenarios(context, ech),
                        get_uncertainties(context, ech),
                        firmness,
                        get_market_schedule(context),
                        get_tso_schedule(context),
                        get_tso_actions(context),
                        gratis_starts,
                        runnable.configs
                        )
end

function update_market_schedule!(context::AbstractContext, ech,
                                result::EnergyMarketModel,
                                firmness,
                                runnable::Union{EnergyMarket,BalanceMarket})
    market_schedule = get_market_schedule(context)
    market_schedule.decider_type = DeciderType(runnable)
    market_schedule.decision_time = ech

    # Pilotables levels
    for ((gen_id, ts, s), p_injected_var) in result.pilotable_model.p_injected
        if get_power_level_firmness(firmness, gen_id, ts) == FREE
            set_prod_value!(market_schedule, gen_id, ts, s, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == TO_DECIDE
            set_prod_definitive_value!(market_schedule, gen_id, ts, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == DECIDED
            @assert( ismissing(get_prod_value(market_schedule, gen_id, ts))
                    || (value(p_injected_var) ≈ get_prod_value(market_schedule, gen_id, ts)) )
            set_prod_definitive_value!(market_schedule, gen_id, ts, value(p_injected_var))
        end
    end

    # Commitment
    for ((gen_id, ts, s), b_on_var) in result.pilotable_model.b_on
        gen_state_value = parse(GeneratorState, value(b_on_var))
        if get_commitment_firmness(firmness, gen_id, ts) == FREE
            set_commitment_value!(market_schedule, gen_id, ts, s, gen_state_value)
        elseif get_commitment_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            set_commitment_definitive_value!(market_schedule, gen_id, ts, gen_state_value)
        end
    end

    # Capping
    update_schedule_capping!(market_schedule, context, ech, result.limitable_model, runnable.configs.CONSIDER_TSOACTIONS_LIMITATIONS)

    # Limitables levels : needs to be after capping update
    uncertainties_l = get_uncertainties(context, ech)
    for gen in Networks.get_generators_of_type(get_network(context), Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        for ts in get_target_timepoints(context)
            for s in get_scenarios(context)
                capped_l = safeget_capping(market_schedule, gen_id, ts, s)
                injected_l = get_uncertainties(uncertainties_l,gen_id,ts,s) - capped_l
                set_prod_value!(market_schedule, gen_id, ts, s, injected_l)
            end
        end
    end

    # loss_of_load (load-shedding)
    update_schedule_loss_of_load!(market_schedule, context, ech, result.lol_model)

    return market_schedule
end

