using .Networks

using Printf
using Dates
using Parameters
using DataStructures

@enum GeneratorState begin
    ON
    OFF
end

function Base.float(generator_state::GeneratorState)
    if generator_state == ON
        return 1.
    elseif generator_state == OFF
        return 0.
    else
        throw( error("Unknown GeneratorState value : ", generator_state) )
    end
end

function Base.parse(::Type{GeneratorState}, str::String)
    if lowercase(str) == "on"
        return ON
    elseif  lowercase(str) == "off"
        return OFF
    else
        throw( error("Unable to convert `", str, "` to a GeneratorState") )
    end
end

function Base.parse(::Type{GeneratorState}, val::Float64)
    if val > 1e-09
        return ON
    else
        return OFF
    end
end

##########################################
## Firmness
##########################################

@enum DecisionFirmness begin
    DECIDED # a firm decision was already taken
    TO_DECIDE # a firm decision needs to be taken
    FREE # by scenario decisions
end

@with_kw_noshow struct Firmness
    #gen,ts
    commitment::SortedDict{String, SortedDict{Dates.DateTime, DecisionFirmness} } =
        SortedDict{String, SortedDict{Dates.DateTime, DecisionFirmness} }()
    power_level::SortedDict{String, SortedDict{Dates.DateTime, DecisionFirmness} } =
        SortedDict{String, SortedDict{Dates.DateTime, DecisionFirmness} }()
end

Base.:(==)(a::Firmness, b::Firmness) = a.commitment==b.commitment && a.power_level==b.power_level

function get_commitment_firmness(firmness::Firmness)
    return firmness.commitment
end

function get_power_level_firmness(firmness::Firmness)
    return firmness.power_level
end

function set_commitment_firmness!(firmness::Firmness, gen_id::String, ts, decision_firmness::DecisionFirmness)
    get!(firmness.commitment, gen_id, SortedDict{Dates.DateTime, DecisionFirmness}())
    firmness.commitment[gen_id][ts] = decision_firmness
end

function set_power_level_firmness!(firmness::Firmness, gen_id::String, ts, decision_firmness::DecisionFirmness)
    get!(firmness.power_level, gen_id, SortedDict{Dates.DateTime, DecisionFirmness}())
    firmness.power_level[gen_id][ts] = decision_firmness
end

function get_commitment_firmness(firmness::Firmness, gen_id::String)
    if haskey(firmness.commitment,gen_id)
        return firmness.commitment[gen_id]
    else
        return missing
    end
end

function get_commitment_firmness(firmness::Firmness, gen_id::String, ts::Dates.DateTime)
    if haskey(firmness.commitment,gen_id)
        if haskey(firmness.commitment[gen_id], ts)
            return firmness.commitment[gen_id][ts]
        else
            return missing
        end
    else
        return missing
    end
end

function get_power_level_firmness(firmness::Firmness, gen_id::String)
    if haskey(firmness.power_level,gen_id)
        return firmness.power_level[gen_id]
    else
        return missing
    end
end

function get_power_level_firmness(firmness::Firmness, gen_id::String, ts::Dates.DateTime)
    if haskey(firmness.power_level,gen_id)
        if haskey(firmness.power_level[gen_id], ts)
            return firmness.power_level[gen_id][ts]
        else
            return missing
        end
    else
        return missing
    end
end

function Base.show(io::IO, firmness::Firmness)
    println(io, "commitment :")
    pretty_print(io, firmness.commitment)
    println(io, "power_level :")
    pretty_print(io, firmness.power_level)
end

##########################################
## Schedule
##########################################

struct GeneratorSchedule
    gen_id::String
    commitment::SortedDict{Dates.DateTime, UncertainValue{GeneratorState}}
    production::SortedDict{Dates.DateTime, UncertainValue{Float64}}
end
function GeneratorSchedule(gen_id::String)
    return  GeneratorSchedule(gen_id,
                            SortedDict{Dates.DateTime, UncertainValue{GeneratorState}}(),
                            SortedDict{Dates.DateTime, UncertainValue{Float64}}())
end


function get_commitment(generator_schedule::GeneratorSchedule)
    return generator_schedule.commitment
end

function get_production(generator_schedule::GeneratorSchedule)
    return generator_schedule.production
end

mutable struct Schedule <: AbstractSchedule
    decider_type::DeciderType
    decision_time::Dates.DateTime
    #gen_id => ts => uncertainValue(s)
    generator_schedules::SortedDict{String, GeneratorSchedule }

    #bus,ts,s
    loss_of_load_by_bus::SortedDict{Tuple{String,DateTime,String}, Float64 }
    #limitable_gen_id,ts,s
    capping::SortedDict{Tuple{String,DateTime,String}, Float64 }
end
function Schedule(decider_type, ech)
    return Schedule(decider_type, ech,
                    SortedDict{String, GeneratorSchedule}(),
                    SortedDict{Tuple{String,DateTime,String}, Float64 }(),
                    SortedDict{Tuple{String,DateTime,String}, Float64 }()
                    )
end
function Schedule(decider_type, ech, generator_schedules)
    return Schedule(decider_type, ech, generator_schedules,
                    SortedDict{Tuple{String,DateTime,String}, Float64 }(),
                    SortedDict{Tuple{String,DateTime,String}, Float64 }()
                    )
end

function init!(schedule::Schedule, network::Networks.Network,
            target_timepoints::Vector{Dates.DateTime}, scenarios::Vector{String})
    empty!(schedule.generator_schedules)
    for generator_l in Networks.get_generators(network)
        gen_id = Networks.get_id(generator_l)
        schedule.generator_schedules[gen_id] = GeneratorSchedule(gen_id)
        for ts in target_timepoints

            schedule.generator_schedules[gen_id].production[ts] = UncertainValue{Float64}(scenarios)

            #commitment values are only defined for generators that have a pmin > 0
            if (Networks.get_p_min(generator_l) > 1e-09)
                schedule.generator_schedules[gen_id].commitment[ts] = UncertainValue{GeneratorState}(scenarios)
            end
        end
    end
end

function get_sub_schedule(schedule::Schedule, gen_id::String)::Union{GeneratorSchedule,Missing}
    if haskey(schedule.generator_schedules, gen_id)
        return schedule.generator_schedules[gen_id]
    else
        throw( error("no generator schedule was initialized for generator ", gen_id))
        #return missing
    end
end

function get_prod_uncertain_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::Union{UncertainValue{Float64},Missing}
    if haskey(sub_schedule.production, ts)
        return sub_schedule.production[ts]
    else
        @warn("power level schedule is not defined for generator ", sub_schedule.gen_id)
        return missing
    end
end
function get_prod_uncertain_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::Union{UncertainValue{Float64},Missing}
    return get_prod_uncertain_value(schedule.generator_schedules[gen_id], ts)
end

function get_commitment_uncertain_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::Union{UncertainValue{GeneratorState},Missing}
    if haskey(sub_schedule.commitment, ts)
        return sub_schedule.commitment[ts]
    else
        # @warn("commitment schedule is not defined for generator ", sub_schedule.gen_id)
        return missing
    end
end
function get_commitment_uncertain_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::Union{UncertainValue{GeneratorState},Missing}
    return get_commitment_uncertain_value(schedule.generator_schedules[gen_id], ts)
end


function get_prod_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::Union{Float64, Missing}
    uncertain_value = get_prod_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value)
end
function get_prod_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::Union{Float64, Missing}
    return get_prod_value(schedule.generator_schedules[gen_id], ts)
end

function safeget_prod_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::Float64
    return safeget_prod_value(schedule.generator_schedules[gen_id], ts)
end
function safeget_prod_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::Float64
    prod = get_prod_value(sub_schedule, ts)
    if ismissing(prod)
        msg = @sprintf("Missing production value for ts=%s in schedule %s",
                        ts, sub_schedule.gen_id)
        throw( error(msg) )
    else
        return prod
    end
end

function get_commitment_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::Union{GeneratorState, Missing}
    uncertain_value = get_commitment_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value)
end
function get_commitment_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::Union{GeneratorState, Missing}
    return get_commitment_value(schedule.generator_schedules[gen_id], ts)
end

function safeget_commitment_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime)::GeneratorState
    return safeget_commitment_value(schedule.generator_schedules[gen_id], ts)
end
function safeget_commitment_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime)::GeneratorState
    commitment = get_commitment_value(sub_schedule, ts)
    if ismissing(commitment)
        msg = @sprintf("Missing commitment value for ts=%s in schedule %s",
                        ts, sub_schedule.gen_id)
        throw( error(msg) )
    else
        return commitment
    end
end

function get_prod_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, scenario::String)::Union{Float64, Missing}
    uncertain_value = get_prod_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value, scenario)
end
function get_prod_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario::String)::Union{Float64, Missing}
    return get_prod_value(schedule.generator_schedules[gen_id], ts, scenario)
end
function safeget_prod_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario::String)::Union{Float64, Missing}
    prod = get_prod_value(schedule.generator_schedules[gen_id], ts, scenario)
    if ismissing(prod)
        msg = @sprintf("Missing production value in schedule for (gen_id=%s,ts=%s,s=%s)",
                        gen_id,ts,scenario)
        throw( error(msg) )
    else
        return prod
    end
end
function safeget_prod_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, scenario::String)::Float64
    prod = get_prod_value(sub_schedule, ts, scenario)
    if ismissing(prod)
        msg = @sprintf("Missing production value for (ts=%s,s=%s) in schedule %s",
                        ts, scenario, sub_schedule.gen_id)
        throw( error(msg) )
    else
        return prod
    end
end

function get_commitment_value(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, scenario::String)::Union{GeneratorState, Missing}
    uncertain_value = get_commitment_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value, scenario)
end
function get_commitment_value(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario::String)::Union{GeneratorState, Missing}
    return get_commitment_value(schedule.generator_schedules[gen_id], ts, scenario)
end

function get_capping(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario)::Union{Float64, Missing}
    key_l = (gen_id, ts, scenario)
    if !haskey(schedule.capping, key_l)
        return missing
    else
        return schedule.capping[key_l]
    end
end
function safeget_capping(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario)::Float64
    capping = get_capping(schedule, gen_id, ts, scenario)
    if ismissing(capping)
        msg = @sprintf("Missing capping value for (gen_id=%s,ts=%s,s=%s) in schedule",
                        gen_id, ts, scenario)
        throw( error(msg) )
    else
        return capping
    end
end

function get_loss_of_load(schedule::Schedule, bus_id, ts::Dates.DateTime, scenario)::Union{Float64, Missing}
    key_l = (bus_id, ts, scenario)
    if !haskey(schedule.loss_of_load_by_bus, key_l)
        return missing
    else
        return schedule.loss_of_load_by_bus[key_l]
    end
end
function safeget_loss_of_load(schedule::Schedule, bus_id, ts::Dates.DateTime, scenario)::Float64
    loss_of_load = get_loss_of_load(schedule, bus_id, ts, scenario)
    if ismissing(loss_of_load)
        msg = @sprintf("Missing loss_of_load value for (gen_id=%s,ts=%s,s=%s) in schedule",
                        bus_id, ts, scenario)
        throw( error(msg) )
    else
        return loss_of_load
    end
end

function set_prod_value!(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, scenario::String, value::Float64)
    uncertain_value = get_prod_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        msg = @sprintf("cannot set uninitialized schedule production value for generator `%s` at timestep %s and scenario %s",
                        sub_schedule.gen_id, ts, scenario)
        throw( error(msg) )
    end
    set_value!(uncertain_value, scenario, value)
end
function set_prod_value!(schedule, gen_id::String, ts::Dates.DateTime, scenario::String, value::Float64)
    return set_prod_value!(schedule.generator_schedules[gen_id], ts, scenario, value)
end

function set_commitment_value!(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, scenario::String, value::GeneratorState)
    uncertain_value = get_commitment_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        msg = @sprintf("cannot set uninitialized schedule commitment value for generator `%s` at timestep %s and scenario %s",
                        sub_schedule.gen_id, ts, scenario)
        throw( error(msg) )
    end
    set_value!(uncertain_value, scenario, value)
end
function set_commitment_value!(schedule, gen_id::String, ts::Dates.DateTime, scenario::String, value::GeneratorState)
    return set_commitment_value!(schedule.generator_schedules[gen_id], ts, scenario, value)
end


function set_prod_definitive_value!(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, value::Float64)
    uncertain_value = get_prod_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        msg = @sprintf("cannot set uninitialized schedule production value for generator `%s` at timestep %s",
                        sub_schedule.gen_id, ts)
        throw( error(msg) )
    end
    set_definitive_value!(uncertain_value, value)
end
function set_prod_definitive_value!(schedule, gen_id::String, ts::Dates.DateTime, value::Float64)
    return set_prod_definitive_value!(schedule.generator_schedules[gen_id], ts, value)
end

function set_commitment_definitive_value!(sub_schedule::GeneratorSchedule, ts::Dates.DateTime, value::GeneratorState)
    uncertain_value = get_commitment_uncertain_value(sub_schedule, ts)
    if ismissing(uncertain_value)
        msg = @sprintf("cannot set uninitialized schedule commitment value for generator `%s` at timestep %s",
                        sub_schedule.gen_id, ts)
        throw( error(msg) )
    end
    set_definitive_value!(uncertain_value, value)
end
function set_commitment_definitive_value!(schedule, gen_id::String, ts::Dates.DateTime, value::GeneratorState)
    return set_commitment_definitive_value!(schedule.generator_schedules[gen_id], ts, value)
end

function set_capping_value!(schedule::Schedule, gen_id::String, ts::Dates.DateTime, scenario::String, value::Float64)
    schedule.capping[gen_id, ts, scenario] = value
end
function set_loss_of_load_value!(schedule::Schedule, bus_id::String, ts::Dates.DateTime, scenario::String, value::Float64)
    schedule.loss_of_load_by_bus[bus_id, ts, scenario] = value
end

function reset_loss_of_load_by_bus!(schedule::Schedule)
    empty!(schedule.loss_of_load_by_bus)
end
function reset_capping!(schedule::Schedule)
    empty!(schedule.capping)
end

function get_commitment_sub_schedule(schedule::Schedule, gen_id::String)
    return get_commitment(get_sub_schedule(schedule, gen_id))
end

function get_production_sub_schedule(schedule::Schedule, gen_id::String)
    return get_production(get_sub_schedule(schedule, gen_id))
end

function Base.show(io::IO, gen_schedule::GeneratorSchedule)
    @printf("generator %s:\n", gen_schedule.gen_id)
    println("commitment:")
    pretty_print(io, gen_schedule.commitment)
    println("production:")
    pretty_print(io, gen_schedule.production)
end
function Base.show(io::IO, schedule::Schedule)
    @printf("schedule decided by %s at %s:\n", schedule.decider_type, schedule.decision_time)
    pretty_print(io, schedule.generator_schedules)
    println("capping:")
    pretty_print(io, schedule.capping)
    println("loss_of_load_by_bus:")
    pretty_print(io, schedule.loss_of_load_by_bus)
end


#######################
# Helpers
########################

function get_target_timepoints(gen_schedule::GeneratorSchedule)
    return union(keys(gen_schedule.commitment), keys(gen_schedule.production))
end

function get_start_value(preceding_state, current_state)
    if preceding_state==OFF && current_state==ON
        return 1
    else
        return 0
    end
end

function get_start_value(generator_reference_schedule, ts, preceding_ts, generator_initial_state)
    current_state = safeget_commitment_value(generator_reference_schedule, ts)
    preceding_state = ( isnothing(preceding_ts) ? generator_initial_state :
                            safeget_commitment_value(generator_reference_schedule, preceding_ts) )
    return get_start_value(preceding_state, current_state)
end
