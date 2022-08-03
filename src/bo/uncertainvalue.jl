using .Networks

using Printf
using DataStructures

mutable struct UncertainValue{T}
    definitive_value::Union{T, Missing}
    anticipated_value::SortedDict{String, Union{T, Missing}}
end

function UncertainValue{T}() where T
    return UncertainValue{T}( missing, SortedDict{String, T}() )
end
function UncertainValue{T}(scenarios::Vector{String}) where T
    antcipated_values = SortedDict(zip(scenarios, fill(missing, length(scenarios))))
    return UncertainValue{T}(missing, antcipated_values)
end

function add_missing_scenarios(uncertain_value::UncertainValue{T}, scenarios) where T
    for scenario in scenarios
        get!(uncertain_value.anticipated_value, scenario, missing)
    end
end

function get_scenarios(uncertain_value::UncertainValue{T}) where T
    return keys(uncertain_value.anticipated_value)
end

function is_missing_values(uncertain_value::UncertainValue{T})::Bool where T
    for (_,value) in uncertain_value.anticipated_value
        if ismissing(value)
            return true
        end
    end
    return false
end

function make_non_definitive(uncertain_value::UncertainValue{T}) where T
    uncertain_value.definitive_value = missing
end

function set_value!(uncertain_value::UncertainValue{T}, scenario::String, value::T)::Union{T, Missing} where T
    if is_definitive(uncertain_value)
        make_non_definitive(uncertain_value)
    end
    uncertain_value.anticipated_value[scenario] = value
end

function is_definitive(uncertain_value::UncertainValue{T})::Bool where T
    existing_value = get_value(uncertain_value)
    return !ismissing(existing_value)
end

function safeset_definitive_value!(uncertain_value::UncertainValue{T}, value::T)::Union{T, Missing} where T
    if is_definitive(uncertain_value)
        existing_value = get_value(uncertain_value)
        if is_different(existing_value, value)
            msg = @sprintf("Unable to set definitive value to %s : The definitive value `%s` was already set to the UncertainValue.",
                            value, existing_value)
            throw( error(msg) )
        end
    else
        uncertain_value.definitive_value = value
        for (scenario, _) in uncertain_value.anticipated_value
            uncertain_value.anticipated_value[scenario] = value
        end
    end
    return uncertain_value.definitive_value
end
function set_definitive_value!(uncertain_value::UncertainValue{T}, value::T)::Union{T, Missing} where T
    if is_definitive(uncertain_value)
        existing_value = get_value(uncertain_value)
        if existing_value != value
            msg = @sprintf("changed value to %s : The definitive value was already set to `%s` for the UncertainValue.",
                    value, existing_value)
            @debug msg
        end
    end

    uncertain_value.definitive_value = value
    for (scenario, _) in uncertain_value.anticipated_value
        uncertain_value.anticipated_value[scenario] = value
    end
    return uncertain_value.definitive_value
end

function get_value(uncertain_value::UncertainValue{T})::Union{T, Missing} where T
    return uncertain_value.definitive_value
end

function safeget_value(uncertain_value::UncertainValue{T})::T where T
    if is_definitive(uncertain_value)
        return get_value(uncertain_value)
    else
        throw( error("UncertainValue is not definitive") )
    end
end

function get_value(uncertain_value::UncertainValue{T}, scenario::String)::Union{Missing,T} where T
    if is_definitive(uncertain_value)
        return get_value(uncertain_value)
    else
        if haskey(uncertain_value.anticipated_value, scenario)
            return uncertain_value.anticipated_value[scenario]
        else
            return missing
        end
    end
end

function safeget_value(uncertain_value::UncertainValue{T}, scenario::String)::T where T
    value = get_value(uncertain_value, scenario)
    if !ismissing( value )
        return value
    else
        msg = @sprintf("No decision was made for scenario %s yet", scenario)
        throw( error(msg) )
    end
end

function Base.show(io::IO, uncertain_value::UncertainValue{T}) where T
    scenarios = get_scenarios(uncertain_value)
    if is_definitive(uncertain_value)
        @printf(io, "definitive value %s for scenarios %s",
                get_value(uncertain_value), collect(scenarios))
    else
        for s in scenarios
            @printf(io, "%s:%s,", s, get_value(uncertain_value, s))
        end
    end
end