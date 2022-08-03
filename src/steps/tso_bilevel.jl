using .Networks

using Dates
using JuMP
using Printf
using Parameters

@with_kw struct TSOBilevel <: AbstractTSO
    configs::TSOBilevelConfigs = TSOBilevelConfigs()
end

function run(runnable::TSOBilevel, ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime}, context::AbstractContext)
    problem_name_l = @sprintf("tso_bilevel_%s", ech)

    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    return tso_bilevel(get_network(context),
                    TS,
                    get_generators_initial_state(context),
                    get_scenarios(context, ech),
                    get_uncertainties(context, ech),
                    firmness,
                    get_market_schedule(context),
                    get_tso_schedule(context),
                    runnable.configs
                    )
end

function update_tso_schedule!(context::AbstractContext, ech, result::TSOBilevelModel, firmness,
                            runnable::TSOBilevel)
    tso_schedule = get_tso_schedule(context)
    tso_schedule.decider_type = DeciderType(runnable)
    tso_schedule.decision_time = ech

    # upper problem (TSO) locates limitable injections
    for ((gen_id, ts, s), p_injected_var) in result.upper.limitable_model.p_injected
        set_prod_value!(tso_schedule, gen_id, ts, s, value(p_injected_var))
    end
    # lower problem (Market) decides pilotable injections
    for ((gen_id, ts, s), p_injected_var) in result.lower.pilotable_model.p_injected
        if get_power_level_firmness(firmness, gen_id, ts) == FREE
            set_prod_value!(tso_schedule, gen_id, ts, s, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == TO_DECIDE
            set_prod_definitive_value!(tso_schedule, gen_id, ts, value(p_injected_var))
        elseif get_power_level_firmness(firmness, gen_id, ts) == DECIDED
            @assert( ismissing(get_prod_value(tso_schedule, gen_id, ts))
                    || value(p_injected_var) ≈ get_prod_value(tso_schedule, gen_id, ts) )
            set_prod_definitive_value!(tso_schedule, gen_id, ts, value(p_injected_var))
        end
    end

    for ((gen_id, ts, s), b_on_var) in result.upper.pilotable_model.b_on
        gen_state_value = parse(GeneratorState, value(b_on_var))
        if get_commitment_firmness(firmness, gen_id, ts) == FREE
            set_commitment_value!(tso_schedule, gen_id, ts, s, gen_state_value)
        elseif get_commitment_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            set_commitment_definitive_value!(tso_schedule, gen_id, ts, gen_state_value)
        end
    end

    # Capping : upper problem (TSO) locates cappings
    update_schedule_capping!(tso_schedule, result.upper.limitable_model)

    # loss_of_load (load-shedding) : upper problem (TSO) locates load shedding
    update_schedule_loss_of_load!(tso_schedule, result.upper.lol_model)

    return tso_schedule
end

function update_schedule_capping!(tso_schedule, limitable_model::TSOBilevelTSOLimitableModel)
    reset_capping!(tso_schedule)

    for ((gen_id, ts, s), p_capping_var) in limitable_model.p_capping
        tso_schedule.capping[gen_id, ts, s] = value(p_capping_var)
    end
end

function update_schedule_loss_of_load!(tso_schedule, lol_model::TSOBilevelTSOLoLModel)
    reset_loss_of_load_by_bus!(tso_schedule)

    for ((bus_id, ts, s), p_loss_of_load_var) in lol_model.p_loss_of_load
        tso_schedule.loss_of_load_by_bus[bus_id, ts, s] = value(p_loss_of_load_var)
    end
end


function update_tso_actions!(context::AbstractContext, ech, result, firmness,
                            runnable::TSOBilevel)
    tso_actions = get_tso_actions(context)
    reset_tso_actions!(tso_actions)

    # Limitations : only firm i.e. value is common to all scenarios
    limitations = SortedDict{Tuple{String,DateTime}, Float64}() #TODELETE
    for ((gen_id, ts, s), p_limit_var) in result.upper.limitable_model.p_limit
        if (value(result.upper.limitable_model.b_is_limited[gen_id, ts, s]) > 1e-09)
            if ( get_power_level_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
                || runnable.configs.LINK_SCENARIOS_LIMIT )
                @assert( value(p_limit_var) ≈ get!(limitations, (gen_id, ts), value(p_limit_var)) ) #TODELETE : checks that all values are the same across scenarios
                add_missing_scenarios(get_limitation_uncertain_value!(tso_actions, gen_id, ts), get_scenarios(context))
                set_limitation_definitive_value!(tso_actions, gen_id, ts, value(p_limit_var))
            else
                #FIXME : may encounter problems if runnable.configs.LINK_SCENARIOS_LIMIT==false, cause limitations are supposed firm
                @warn "FIXME? : need to fix limitation actions to a by scenario before DP"
                set_limitation_value!(tso_actions, gen_id, ts, s, value(p_limit_var))
            end
        #else : will remain missing
        end
    end

    # Impositions
    impositions = SortedDict{Tuple{String,DateTime}, Float64}() #TODELETE
    for ((gen_id, ts, s), p_injected_var) in result.lower.pilotable_model.p_injected
        p_min_var = result.upper.pilotable_model.p_imposition_min[gen_id, ts, s]
        p_max_var = result.upper.pilotable_model.p_imposition_max[gen_id, ts, s]

        if get_power_level_firmness(firmness, gen_id, ts) in [TO_DECIDE, DECIDED]
            @assert( value(p_injected_var) ≈ get!(impositions, (gen_id, ts), value(p_injected_var)) ) #TODELETE : checks that all values are the same across scenarios
            set_imposition_value!(tso_actions, gen_id, ts, s, value(p_min_var), value(p_max_var))
        else
            set_imposition_value!(tso_actions, gen_id, ts, s, value(p_min_var), value(p_max_var))

        end
    end

end
