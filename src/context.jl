using Dates

using ..Networks

mutable struct PSCOPFContext <: AbstractContext
    network::Networks.Network
    target_timepoints::Vector{Dates.DateTime}
    management_mode::ManagementMode

    #FIXME : Question : Do we need an initial schedule for power levels ?
    #Not if we make sure ECH[1] <= TS-DMO

    generators_initial_state::SortedDict{String,GeneratorState}

    #uncertainties
    uncertainties::Uncertainties

    #AssessmentUncertainties
    assessment_uncertainties

    horizon_timepoints::Vector{Dates.DateTime}

    market_schedule::Schedule
    tso_schedule::Schedule

    tso_actions::TSOActions
    #flows ?

    out_dir
end

function PSCOPFContext(network::Networks.Network, target_timepoints::Vector{Dates.DateTime},
                    management_mode::ManagementMode,
                    generators_initial_state::SortedDict{String,GeneratorState}=SortedDict{String,GeneratorState}(),
                    uncertainties::Uncertainties=Uncertainties(),
                    assessment_uncertainties=nothing,
                    out_dir=nothing
                    )
    market_schedule = Schedule(Utilitary(), Dates.DateTime(0))
    init!(market_schedule, network, target_timepoints, get_scenarios(uncertainties))
    tso_schedule = Schedule(Utilitary(), Dates.DateTime(0))
    init!(tso_schedule, network, target_timepoints, get_scenarios(uncertainties))
    return PSCOPFContext(network, target_timepoints, management_mode,
                        generators_initial_state,
                        uncertainties, assessment_uncertainties,
                        Vector{Dates.DateTime}(),
                        market_schedule,
                        tso_schedule,
                        TSOActions(),
                        out_dir)
end

function get_network(context::PSCOPFContext)
    return context.network
end

function get_target_timepoints(context::PSCOPFContext)
    return context.target_timepoints
end

function get_horizon_timepoints(context::PSCOPFContext)
    return context.horizon_timepoints
end

function set_horizon_timepoints(context::PSCOPFContext, horizon_timepoints::Vector{Dates.DateTime})
    context.horizon_timepoints = horizon_timepoints
end

function get_management_mode(context::PSCOPFContext)
    return context.management_mode
end

function get_generators_initial_state(context::PSCOPFContext)
    return context.generators_initial_state
end

function get_uncertainties(context::PSCOPFContext)::Uncertainties
    return context.uncertainties
end

function get_uncertainties(context::PSCOPFContext, ech::Dates.DateTime)::UncertaintiesAtEch
    return get_uncertainties(context)[ech]
end

function get_scenarios(context::PSCOPFContext, ech::Dates.DateTime)::Vector{String}
    uncertainties_at_ech = get_uncertainties(context, ech)
    if isnothing(uncertainties_at_ech)
        return Vector{String}()
    else
        return get_scenarios(uncertainties_at_ech)
    end
end

function get_scenarios(context::PSCOPFContext)::Vector{String}
    uncertainties = get_uncertainties(context)
    if isnothing(uncertainties)
        return Vector{String}()
    else
        return get_scenarios(uncertainties)
    end
end

function get_assessment_uncertainties(context::PSCOPFContext)
    return context.assessment_uncertainties
end

function set_current_ech!(context_p::PSCOPFContext, ech::Dates.DateTime)
    context_p.current_ech = ech
end

function get_tso_schedule(context_p::PSCOPFContext)
    return context_p.tso_schedule
end

function get_market_schedule(context_p::PSCOPFContext)
    return context_p.market_schedule
end


function get_limitables_ids(context_p::PSCOPFContext)
    limitables = Networks.get_generators_of_type(get_network(context_p), Networks.LIMITABLE)
    limitables_ids = map(lim_gen->Networks.get_id(lim_gen), limitables)
    return limitables_ids
end

function get_initial_state(initial_states::SortedDict{String,GeneratorState}, generator::Generator)
    if Networks.get_p_min(generator) < 1e-09
        return ON
    else
        return initial_states[Networks.get_id(generator)]
    end
end

function definitive_starts(schedule::Schedule, initial_state::SortedDict{String, GeneratorState})
    result = Set{Tuple{String,Dates.DateTime}}()

    for (gen_id, gen_schedule) in schedule.generator_schedules
        if isempty(gen_schedule.commitment)
            continue
        end

        prev_ts = nothing
        prev_state = initial_state[gen_id]

        for (ts, current_state) in gen_schedule.commitment
            @assert( isnothing(prev_ts) || (prev_ts < ts) )
            if !is_definitive(current_state)
                break #if current state is not definitive the following are not neither

            elseif ( (prev_state==OFF) && (get_value(current_state)==ON) )
                push!(result, (gen_id,ts) )
                prev_state = get_value(current_state)
                prev_ts = ts
            end
        end
    end

    return result
end
