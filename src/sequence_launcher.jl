using Dates
using DataStructures
using Parameters
using Printf

using .Networks

###############################################
# Sequence
###############################################

"""
    Sequence

Lists, for each horizon timepoint, the ordered operations to execute at that time.
"""
@with_kw_noshow struct Sequence
    operations::SortedDict{Dates.DateTime, Vector{AbstractRunnable}} = SortedDict{Dates.DateTime, Vector{AbstractRunnable}}()
end

function get_operations(sequence::Sequence)
    return sequence.operations
end

function get_timepoints(sequence::Sequence)
    return collect(keys(get_operations(sequence)))
end

function get_steps(sequence::Sequence, ech::DateTime)
    return get_operations(sequence)[ech]
end

"""
    length of the sequence in terms of time not in terms of number of runnables to execute
"""
function Base.length(sequence::Sequence)
    return length(get_operations(sequence))
end

function Base.show(io::IO, sequence::Sequence)
    pretty_print(io, get_operations(sequence))
end

function get_horizon_timepoints(sequence::Sequence)::Vector{Dates.DateTime}
    return collect(keys(get_operations(sequence)))
end

function get_ech(sequence::Sequence, index::Int)
    if length(sequence) < index
        throw( error("attempt to acess ", length(sequence), "-element Sequence at index ", index, ".") )
    else
        return get_horizon_timepoints(sequence)[index]
    end
end

function add_step!(sequence::Sequence, step_type::Type{T}, ech::Dates.DateTime) where T<:AbstractRunnable
    step_instance = step_type()
    add_step!(sequence, step_instance, ech)
end
function add_step!(sequence::Sequence, step_instance::AbstractRunnable, ech::Dates.DateTime)
    steps_at_ech = get!(sequence.operations, ech, Vector{AbstractRunnable}())
    push!(steps_at_ech, step_instance)
end

###############################################
# Sequence Launch
###############################################

function init!(context_p::AbstractContext, sequence_p::Sequence, check_context::Bool)
    if !isnothing(context_p.out_dir)
        @warn("removing files that start with $(PSCOPFio.OUTPUT_PREFIX) in $(context_p.out_dir)")
        rm_prefixed(context_p.out_dir, PSCOPFio.OUTPUT_PREFIX)
    end

    println("Lancement du mode : ", context_p.management_mode.name)
    println("Dates d'interet : ", get_target_timepoints(context_p))
    set_horizon_timepoints(context_p, get_timepoints(sequence_p))
    println("Dates d'échéances : ", get_horizon_timepoints(context_p))

    if check_context && !check(context_p)
        throw( error("Invalid context!") )
    end
end

function run_step!(context_p::AbstractContext, step::AbstractRunnable, ech, next_ech)
    println("-"^20)
    println(typeof(step), " à l'échéance ", ech)
    firmness = compute_firmness(step, ech, next_ech,
                            get_target_timepoints(context_p), context_p)
    trace_firmness(firmness)
    result = run(step, ech, firmness,
                get_target_timepoints(context_p),
                context_p)

    if affects_market_schedule(step)
        println("update market schedule based on optimization results")
        old_market_schedule = deepcopy(get_market_schedule(context_p))
        update_market_schedule!(context_p, ech, result, firmness, step)
        println("Changes to the market schedule:")
        trace_delta_schedules(old_market_schedule, get_market_schedule(context_p))
        #TODO : error if !verify
        verify_firmness(firmness, context_p.market_schedule,
                        excluded_ids=get_limitables_ids(context_p))
        PSCOPF.PSCOPFio.write(context_p, get_market_schedule(context_p), "market_")
        PSCOPF.PSCOPFio.write(context_p, get_market_schedule(context_p), get_uncertainties(context_p), "market_")
        update_market_flows!(context_p)
        trace_flows(get_market_flows(context_p), get_network(context_p))
    end

    if affects_tso_schedule(step)
        println("update TSO schedule based on optimization results")
        old_tso_schedule = deepcopy(get_tso_schedule(context_p))
        update_tso_schedule!(context_p, ech, result, firmness, step)
        println("Changes to the TSO schedule:")
        trace_delta_schedules(old_tso_schedule, get_tso_schedule(context_p))
        #TODO : error if !verify
        verify_firmness(firmness, context_p.tso_schedule,
                        excluded_ids=get_limitables_ids(context_p))
        PSCOPF.PSCOPFio.write(context_p, get_tso_schedule(context_p), "tso_")
        PSCOPF.PSCOPFio.write(context_p, get_tso_schedule(context_p), get_uncertainties(context_p), "tso_")
        update_tso_flows!(context_p)
        trace_flows(get_tso_flows(context_p), get_network(context_p))
    end

    if affects_tso_actions(step)
        println("update TSO actions based on optimization results")
        update_tso_actions!(context_p,
                            ech, result, firmness, step)
        trace_tso_actions(get_tso_actions(context_p))
    end
    #TODO check coherence between tso schedule and actions

    if (affects_market_schedule(step) || affects_tso_schedule(step))
        println("Changes between steps:")
        schedules_to_delta = sort([get_market_schedule(context_p), get_tso_schedule(context_p)],
                                    by=x->x.decision_time)
        trace_delta_schedules(schedules_to_delta...)
    end
end

function run!(context_p::AbstractContext, sequence_p::Sequence;
                check_context=true)
    println(sequence_p)
    init!(context_p, sequence_p, check_context)

    for (steps_index, (ech, steps_at_ech)) in enumerate(get_operations(sequence_p))
        println("-"^50)
        delta = Dates.value(Dates.Minute(get_target_timepoints(context_p)[1]-ech))
        println("ECH : ", ech, " : M-", delta)
        println("-"^50)
        for step in steps_at_ech
            next_ech = get_next_ech(sequence_p, steps_index, step)
            run_step!(context_p, step, ech, next_ech)
        end
    end
end


"""
returns the ech that come after ech at which we execute a step of the same DeciderType as the input step
if none exists returns nothing
"""
function get_next_ech(sequence::Sequence, index::Int, decider_step::AbstractRunnable)
    if index >= length(sequence)
        return nothing
    end

    for ech_l in get_horizon_timepoints(sequence)[index+1:end]
        for step_l in get_steps(sequence, ech_l)
            if DeciderType(step_l) == DeciderType(decider_step)
                return ech_l
            end
        end
    end

    return nothing
end


##################################
# Trace
##################################

function trace_firmness(firmness::Firmness)
    for (gen_id,_) in firmness.commitment
        TS_to_decide_for = Vector{DateTime}()
        for (ts, decision_firmness) in firmness.commitment[gen_id]
            if decision_firmness == TO_DECIDE
                push!(TS_to_decide_for, ts)
            end
        end
        if !isempty(TS_to_decide_for)
            @printf("Firm commitment decision must be issued for generator %s for timesteps %s\n",
                    gen_id, TS_to_decide_for)
        end
    end

    for (gen_id,_) in firmness.power_level
        TS_to_decide_for = Vector{DateTime}()
        for (ts, decision_firmness) in firmness.power_level[gen_id]
            if decision_firmness == TO_DECIDE
                push!(TS_to_decide_for, ts)
            end
        end
        if !isempty(TS_to_decide_for)
            @printf("Firm production level decision must be issued for generator %s at timesteps %s\n",
                    gen_id, TS_to_decide_for)
        end
    end
end

function trace_flows(flows::SortedDict{Tuple{String, DateTime, String}, Float64},
                    network::Networks.Network)
    for ((branch_id, ts, scenario), flow_val) in flows
        branch = Networks.safeget_branch(network, branch_id)
        limit = Networks.get_limit(branch)
        if abs(flow_val) > limit
            @printf("Flow value %f for branch %s, at timestep %s and scenario %s exceeds branch limit (%f)\n",
                    flow_val, branch_id, ts, scenario, limit)
        end
    end
end

function trace_delta_schedules(old_schedule::Schedule, new_schedule::Schedule)
    print_non_firm_changes = false

    @printf("Commitment updates : new %s Vs %s\n",
            new_schedule.decider_type, old_schedule.decider_type)
    trace_delta_schedule_component(old_schedule, new_schedule, get_commitment_sub_schedule, print_non_firm_changes)
    @printf("Production updates : new %s Vs %s\n",
            new_schedule.decider_type, old_schedule.decider_type)
    trace_delta_schedule_component(old_schedule, new_schedule, get_production_sub_schedule, print_non_firm_changes)
end

function trace_delta_schedule_component(old_schedule::Schedule, new_schedule::Schedule,
                                        component_accessor::Function,
                                        print_non_firm_changes=false)
    if !print_non_firm_changes
        println("(only showing non firm changes)")
    end
    for (gen_id, _) in new_schedule.generator_schedules
        @assert(haskey(old_schedule.generator_schedules, gen_id))
        old_schedule_component = component_accessor(old_schedule, gen_id)
        new_schedule_component = component_accessor(new_schedule, gen_id)

        for (ts, new_uncertain) in new_schedule_component
            @assert(haskey(old_schedule_component, ts))
            old_uncertain = old_schedule_component[ts]
            if is_definitive(new_uncertain) && is_definitive(old_uncertain)
                if is_different(get_value(new_uncertain), get_value(old_uncertain))
                    @printf("%s changed %s old firm value for generator %s at timestep %s from %s to %s\n",
                            new_schedule.decider_type, old_schedule.decider_type,
                            gen_id, ts, get_value(old_uncertain), get_value(new_uncertain))
                end
            elseif is_definitive(new_uncertain) && !is_definitive(old_uncertain)
                @printf("%s set a firm value for generator %s at timestep %s to %s\n",
                        new_schedule.decider_type, gen_id, ts, get_value(new_uncertain))
            elseif !is_definitive(new_uncertain) && is_definitive(old_uncertain)
                @printf("%s changed the %s firm value for generator %s at timestep %s \
                        from %s to a by scenario value:\n",
                        new_schedule.decider_type, old_schedule.decider_type, gen_id, ts, get_value(new_uncertain))
                println(new_uncertain.anticipated_value)
            elseif print_non_firm_changes
                for scenario in get_scenarios(new_uncertain)
                    old_val = get_value(old_uncertain, scenario)
                    new_val = get_value(new_uncertain, scenario)
                    if ( (ismissing(new_val) != ismissing(old_val))
                        || ( !ismissing(old_val) && !ismissing(new_val) && is_different(new_val, old_val) )
                        )
                        @printf("gen_id=%s, ts=%s, scenario=%s \
                                old_value=%s, new_value=%s\n",
                                gen_id, ts, scenario, old_val, new_val)
                    end
                end
            end
        end
    end
end

function trace_tso_actions(tso_actions::TSOActions)
    println("TSO limitation actions :")
    trace_limitations(tso_actions)
    println("TSO imposition actions :")
    trace_impositions(tso_actions)
end
function trace_limitations(tso_actions)
    for ((gen_id,ts), val_l) in get_limitations(tso_actions)
        println("\tlimited generator %s for timestep %s to %s", gen_id, ts, val_l)
    end
end
function trace_impositions(tso_actions)
    for ((gen_id,ts), (val_min_l,val_max_l)) in get_impositions(tso_actions)
        println("\timposed generator %s for timestep %s to [%s,%s]", gen_id, ts, val_min_l,val_max_l)
    end
end