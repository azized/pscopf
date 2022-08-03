using Dates
using Parameters
using DataStructures

##########################################
## TSOActions
##########################################

Limitations = SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Float64}}
Impositions = SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Tuple{Float64,Float64}}}

@with_kw struct TSOActions
    #gen_id, ts => uncertain Float value
    # Max allowed P of a limitable for a given ts
    limitations::Limitations = Limitations()
    # Imposed P bounds of an pilotable for a given ts
    impositions::Impositions = Impositions()
end

function reset_tso_actions!(tso_actions::TSOActions)
    empty!(tso_actions.limitations)
    empty!(tso_actions.impositions)

    return tso_actions
end

## Limitations
#--------------

function get_limitations(tso_actions::TSOActions)
    return tso_actions.limitations
end


function get_limitation_uncertain_value(tso_actions::TSOActions,
                                        gen_id::String, ts::Dates.DateTime
                                        )::Union{UncertainValue{Float64},Missing}
    limitations = get_limitations(tso_actions)
    return get_limitation_uncertain_value(limitations, gen_id, ts)
end
function get_limitation_uncertain_value(limitations::SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Float64}},
                                        gen_id::String, ts::Dates.DateTime
                                        )::Union{UncertainValue{Float64},Missing}
    if haskey(limitations, (gen_id, ts))
        return limitations[gen_id, ts]
    else
        # msg = @sprintf("limitations is not defined for (%s,%s)", gen_id, ts)
        # @warn msg
        return missing
    end
end

function get_limitation_uncertain_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Union{UncertainValue{Float64},Missing}
    limitations = get_limitations(tso_actions)
    if haskey(limitations, (gen_id, ts))
        return limitations[gen_id, ts]
    else
        limitations[gen_id, ts] = UncertainValue{Float64}()
        return limitations[gen_id, ts]
    end
end

function get_limitation(limitations::Limitations,
                        gen_id::String, ts::Dates.DateTime
                        )::Union{Float64,Missing}
    uncertain_value = get_limitation_uncertain_value(limitations, gen_id, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value)
end
function get_limitation(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Union{Float64,Missing}
    return get_limitation(get_limitations(tso_actions), gen_id, ts)
end

function safeget_limitation(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Float64
    limitation = get_limitation(tso_actions, gen_id, ts)
    if ismissing(limitation)
        msg = @sprintf("Missing limitation value for gen_id=%s,ts=%s",
                        gen_id, ts)
        throw( error(msg) )
    else
        return limitation
    end
end


function get_limitation(limitations::SortedDict{Tuple{String, Dates.DateTime}, UncertainValue{Float64}},
                        gen_id::String, ts::Dates.DateTime, scenario::String
                        )::Union{Float64, Missing}
    uncertain_value = get_limitation_uncertain_value(limitations, gen_id, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value, scenario)
end
function get_limitation(tso_actions::TSOActions,
                        gen_id::String, ts::Dates.DateTime, scenario::String
                        )::Union{Float64, Missing}
    limitations = get_limitations(tso_actions)
    return get_limitation(limitations, gen_id, ts, scenario)
end

function safeget_limitation(tso_actions, gen_id::String, ts::Dates.DateTime, scenario::String)::Float64
    limitation = get_limitation(tso_actions, gen_id, ts, scenario)
    if ismissing(limitation)
        msg = @sprintf("No limitation entry in TSOActions for (gen_id=%s,ts=%s,s=%s).", gen_id, ts, scenario)
        throw(error(msg))
    else
        return limitation
    end
end

function set_limitation_definitive_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, value::Float64)
    uncertain_value = get_limitation_uncertain_value!(tso_actions, gen_id, ts)
    set_definitive_value!(uncertain_value, value)
end

function set_limitation_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, scenario::String, value::Float64)
    uncertain_value = get_limitation_uncertain_value!(tso_actions, gen_id, ts)
    set_value!(uncertain_value, scenario, value)
end

## Impositions
#--------------

function get_impositions(tso_actions::TSOActions)
    return tso_actions.impositions
end

function get_imposition_uncertain_value(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Union{UncertainValue{Tuple{Float64,Float64}},Missing}
    impositions = get_impositions(tso_actions)
    if haskey(impositions, (gen_id, ts))
        return impositions[gen_id, ts]
    else
        # msg = @sprintf("impositions is not defined for (%s,%s)", gen_id, ts)
        # @warn msg
        return missing
    end
end

function get_imposition_uncertain_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Union{UncertainValue{Tuple{Float64,Float64}},Missing}
    impositions = get_impositions(tso_actions)
    if haskey(impositions, (gen_id, ts))
        return impositions[gen_id, ts]
    else
        impositions[gen_id, ts] = UncertainValue{Tuple{Float64,Float64}}()
        return impositions[gen_id, ts]
    end
end

function get_imposition(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Union{Tuple{Float64,Float64},Missing}
    uncertain_value = get_imposition_uncertain_value(tso_actions, gen_id, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value)
end

function safeget_imposition(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime)::Tuple{Float64,Float64}
    imposition = get_imposition(tso_actions, gen_id, ts)
    if ismissing(imposition)
        msg = @sprintf("Missing imposition value for gen_id=%s,ts=%s",
                        gen_id, ts)
        throw( error(msg) )
    else
        return imposition
    end
end

function get_imposition(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, scenario::String)::Union{Tuple{Float64,Float64}, Missing}
    uncertain_value = get_imposition_uncertain_value(tso_actions, gen_id, ts)
    if ismissing(uncertain_value)
        return missing
    end
    return get_value(uncertain_value, scenario)
end

function safeget_imposition(tso_actions, gen_id::String, ts::Dates.DateTime, scenario::String)::Tuple{Float64,Float64}
    imposition = get_imposition(tso_actions, gen_id, ts, scenario)
    if ismissing(imposition)
        msg = @sprintf("No imposition entry in TSOActions for (gen_id=%s,ts=%s,s=%s).", gen_id, ts, scenario)
        throw(error(msg))
    else
        return imposition
    end
end

function get_imposition_level(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, scenario::String)::Union{Float64, Missing}
    imposition = get_imposition(tso_actions, gen_id, ts, scenario)
    if ismissing(imposition)
        return missing
    else
        value_min, value_max = imposition
        if value_min != value_max
            msg = @sprintf("TSOActions for (gen_id=%s,ts=%s,s=%s) imposes interval [%s,%s] and not a single power level.",
                            gen_id, ts, scenario, value_min, value_max)
            throw(error(msg))
        else
            return value_min
        end
    end
end

function safeget_imposition_level(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, scenario::String)::Float64
    imposition = get_imposition_level(tso_actions, gen_id, ts, scenario)
    if ismissing(imposition)
        msg = @sprintf("TSOActions has no imposition value for (gen_id=%s,ts=%s,s=%s).", gen_id, ts, scenario)
        throw(error(msg))
    else
        return imposition
    end
end

function set_imposition_definitive_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, value_min::Float64, value_max::Float64)
    uncertain_value = get_imposition_uncertain_value!(tso_actions, gen_id, ts)
    set_definitive_value!(uncertain_value, (value_min, value_max))
end

function set_imposition_value!(tso_actions::TSOActions, gen_id::String, ts::Dates.DateTime, scenario::String, value_min::Float64, value_max::Float64)
    uncertain_value = get_imposition_uncertain_value!(tso_actions, gen_id, ts)
    set_value!(uncertain_value, scenario, (value_min, value_max))
end

#######################
# Helpers
########################

"""
Performs a partial shallow copy of the input TSOActions
Modifying the returned tso_actions' kept attributes will modify the original one's.
"""
function filter_tso_actions(tso_actions::TSOActions;
                        keep_limitations::Bool=false,
                        keep_impositions::Bool=false)::TSOActions
    limitations_l = keep_limitations ? tso_actions.limitations : SortedDict{Tuple{String, Dates.DateTime}, Float64}()
    impositions_l = keep_impositions ? tso_actions.impositions : SortedDict{Tuple{String, Dates.DateTime}, Tuple{Float64,Float64}}()

    return TSOActions(limitations_l, impositions_l)
end
