using .Networks

using Dates
using JuMP
using Printf
using Parameters

@with_kw struct TSOOutFO <: AbstractTSO
    configs::TSOConfigs = TSOConfigs()
end

TSOInFO = TSOOutFO

function run(runnable::TSOOutFO, ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime}, context::AbstractContext)
    problem_name_l = @sprintf("tso_out_fo_%s", ech)

    tso_starts = get_starts(get_tso_schedule(context), get_generators_initial_state(context))
    market_starts = get_starts(get_market_schedule(context), get_generators_initial_state(context))
    gratis_starts = union(tso_starts, market_starts)

    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    return tso_out_fo(get_network(context),
                    TS,
                    get_generators_initial_state(context),
                    get_scenarios(context, ech),
                    get_uncertainties(context, ech),
                    firmness,
                    get_market_schedule(context),
                    get_tso_schedule(context),
                    gratis_starts,
                    runnable.configs
                    )
end

function update_tso_schedule!(context::AbstractContext, ech, result, firmness,
                            runnable::TSOOutFO)
    tso_schedule = get_tso_schedule(context)
    tso_schedule.decider_type = DeciderType(runnable)
    tso_schedule.decision_time = ech

    for ((gen_id, ts, s), p_injected_var) in result.limitable_model.p_injected
        set_prod_value!(tso_schedule, gen_id, ts, s, value(p_injected_var))
    end
    for ((gen_id, ts, s), p_injected_var) in result.pilotable_model.p_injected
        if get_power_level_firmness(firmness, gen_id, ts) == FREE
            set_prod_value!(tso_schedule, gen_id, ts, s, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == TO_DECIDE
            set_prod_definitive_value!(tso_schedule, gen_id, ts, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == DECIDED
            # NOTE : decided reference may not be the TSO schedule
            @assert( ismissing(get_prod_value(tso_schedule, gen_id, ts))
                    || value(p_injected_var) ≈ get_prod_value(tso_schedule, gen_id, ts) )
            set_prod_definitive_value!(tso_schedule, gen_id, ts, value(p_injected_var))
        end
    end

    for ((gen_id, ts, s), b_on_var) in result.pilotable_model.b_on
        gen_state_value = parse(GeneratorState, value(b_on_var))
        if get_commitment_firmness(firmness, gen_id, ts) == FREE
            set_commitment_value!(tso_schedule, gen_id, ts, s, gen_state_value)
        elseif get_commitment_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            set_commitment_definitive_value!(tso_schedule, gen_id, ts, gen_state_value)
        end
    end

    # Capping
    update_schedule_capping!(tso_schedule, context, ech, result.limitable_model)

    # loss_of_load (load-shedding)
    update_schedule_loss_of_load!(tso_schedule, context, ech, result.lol_model)

    return tso_schedule
end

function update_schedule_capping!(tso_schedule, context, ech,
                                    limitable_model::TSOLimitableModel)
    reset_capping!(tso_schedule)
    for ((gen_id,ts, s), p_injected_var) in limitable_model.p_injected
        available_prod = get_uncertainties(get_uncertainties(context, ech), gen_id, ts, s)
        injected_prod = value(p_injected_var)
        tso_schedule.capping[gen_id, ts, s] = available_prod - injected_prod
    end
end

function update_schedule_loss_of_load!(tso_schedule, context, ech, lol_model::TSOLoLModel)
    reset_loss_of_load_by_bus!(tso_schedule)

    for ((bus_id, ts, s), p_loss_of_load_var) in lol_model.p_loss_of_load
        tso_schedule.loss_of_load_by_bus[bus_id, ts, s] = value(p_loss_of_load_var)
    end
end


function update_tso_actions!(context::AbstractContext, ech, result, firmness,
                            ::TSOOutFO)
    tso_actions = get_tso_actions(context)
    reset_tso_actions!(tso_actions)

    # Limitations :
    # FIXME ; limit only if there is a limitation needed !?
    limitations = SortedDict{Tuple{String,DateTime}, Float64}() #TODELETE
    for ((gen_id, ts, s), p_limit_var) in result.limitable_model.p_limit
        if get_power_level_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            @assert( value(p_limit_var) ≈ get!(limitations, (gen_id, ts), value(p_limit_var)) ) #TODELETE : checks that all values are the same across scenarios
            add_missing_scenarios(get_limitation_uncertain_value!(tso_actions, gen_id, ts), get_scenarios(context))
            set_limitation_definitive_value!(tso_actions, gen_id, ts, value(p_limit_var))
        end
    end

    # Impositions
    impositions = SortedDict{Tuple{String,DateTime}, Float64}() #TODELETE
    for ((gen_id, ts, s), p_injected_var) in result.pilotable_model.p_injected
        if get_power_level_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            @assert( value(p_injected_var) ≈ get!(impositions, (gen_id, ts), value(p_injected_var)) ) #TODELETE : checks that all values are the same across scenarios
            set_imposition_value!(tso_actions, gen_id, ts, s, value(p_injected_var), value(p_injected_var))
        end
    end

end
