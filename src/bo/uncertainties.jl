using Dates
using DataStructures

"""
    UncertaintyErrorNDistribution

describes the uncertainty distribution of an injection (load or renewable)
We consider `ts` as the target timepoint at which the realisation will have place,
and `ech` is the observation point.
The farthest we are from ts, the greater is the uncertainty: we assume an increase  by a ratio of `time_factor` for each hour.
i.e.: factor_l = max{ 0, ( ts - ech )/1h }
adjusted_sigma = factor_l * error_sigma
random_injection = mu * (1+ Normal(0.,  adjusted_sigma))

# Arguments
    - `min_value::Float64` : minimum value that can be assigned to the injection
    - `max_value::Float64` : maxiimum value that can be assigned to the injection
    - `mu::Float64` : average injection value
    - `sigma::Float64` : sigma of the error (0.1 => <=10% of error around mu 70% of the time)
"""
struct UncertaintyErrorNDistribution
    id::String
    min_value::Float64
    max_value::Float64
    mu::Float64
    error_sigma::Float64
end


"""
    UncertaintyNDistribution

describes the uncertainty distribution of an injection (load or renewable)
We consider `ts` as the target timepoint at which the realisation will have place,
and `ech` is the observation point.
The farthest we are from ts, the greater is the uncertainty: we assume an increase  by a ratio of `time_factor` for each hour.
i.e.: factor_l = max{ 0, ( ts - ech )/1h } * time_factor
adjusted_sigma = factor_l * sigma
random_injection = Normal( mu,  adjusted_sigma)

# Arguments
    - `min_value::Float64` : minimum value that can be assigned to the injection
    - `max_value::Float64` : maxiimum value that can be assigned to the injection
    - `mu::Float64` : average injection value
    - `sigma::Float64` : base sigma value to be used as a standard deviation value for the injection uncertainty distribution
    - `time_factor::Float64` : per hour increase ratio of the uncertainty's standard deviation
"""
struct UncertaintyNDistribution
    id::String
    min_value::Float64
    max_value::Float64
    mu::Float64
    sigma::Float64
    time_factor::Float64
end

UncertaintiesDistribution = SortedDict{String, UncertaintyNDistribution}

function add_uncertainty_distribution!(uncertainties_distribution::UncertaintiesDistribution,
                        id::String,
                        min_value::Float64, max_value::Float64,
                        mu::Float64, sigma::Float64,
                        time_factor::Float64=1.)
    uncertainties_distribution[id] = UncertaintyNDistribution(id, min_value, max_value,
                                                            mu, sigma,
                                                            time_factor)
end


# uncertainties[ech][nodal_injection_name][ts][scenario_name]
InjectionUncertainties = SortedDict{Dates.DateTime,
                            SortedDict{String, Float64}
                            }
UncertaintiesAtEch = SortedDict{String, InjectionUncertainties}
Uncertainties = SortedDict{Dates.DateTime, UncertaintiesAtEch}


function add_uncertainty!(uncertainties_at_ech::UncertaintiesAtEch, nodal_injection_name::String, ts::Dates.DateTime, scenario_name::String, val::Float64)
    get!(uncertainties_at_ech, nodal_injection_name, SortedDict{Dates.DateTime, SortedDict{String, Float64}}())
    get!(uncertainties_at_ech[nodal_injection_name], ts, SortedDict{String, Float64}())

    uncertainties_at_ech[nodal_injection_name][ts][scenario_name] = val
    return uncertainties_at_ech
end
function add_uncertainty!(uncertainties::Uncertainties, ech::Dates.DateTime, nodal_injection_name::String, ts::Dates.DateTime, scenario_name::String, val::Float64)
    uncertainties_at_ech = get!(uncertainties, ech, SortedDict{String, SortedDict{Dates.DateTime, SortedDict{String, Float64}}}())
    add_uncertainty!(uncertainties_at_ech, nodal_injection_name, ts, scenario_name, val)
    return uncertainties
end

function get_scenarios(uncertainties_at_ech::UncertaintiesAtEch)::Vector{String}
    for (_, uncertainties_by_ts) in uncertainties_at_ech
        for (_, by_scenario) in uncertainties_by_ts
            return collect(keys(by_scenario))
        end
    end
    return Vector{String}()
end

function get_scenarios(uncertainties::Uncertainties)::Vector{String}
    for (_, uncertainties_at_ech) in uncertainties
        return get_scenarios(uncertainties_at_ech)
    end
    return Vector{String}()
end

function get_target_timepoints(uncertainties_at_ech::UncertaintiesAtEch)::Vector{Dates.DateTime}
    for (_, by_ts) in uncertainties_at_ech
        return collect(keys(by_ts))
    end
    return Vector{Dates.DateTime}()
end
function get_target_timepoints(uncertainties::Uncertainties)::Vector{Dates.DateTime}
    for (_, uncertainties_at_ech) in uncertainties
        return get_target_timepoints(uncertainties_at_ech)
    end
    return Vector{Dates.DateTime}()
end

function get_horizon_timepoints(uncertainties::Uncertainties)::Vector{Dates.DateTime}
    return collect(keys(uncertainties))
end

function get_uncertainties(uncertainties::Uncertainties, ech::Dates.DateTime)::UncertaintiesAtEch
    return uncertainties[ech]
end

function get_uncertainties(uncertainties_at_ech::UncertaintiesAtEch, injection_name::String)::InjectionUncertainties
    return uncertainties_at_ech[injection_name]
end
function get_uncertainties(uncertainties::Uncertainties, ech::Dates.DateTime, injection_name::String)::InjectionUncertainties
    uncertainties_at_ech = get_uncertainties(uncertainties, ech)
    return get_uncertainties(uncertainties_at_ech, injection_name)
end


function get_uncertainties(injection_uncertainties::InjectionUncertainties, ts::Dates.DateTime)
    return injection_uncertainties[ts]
end
function get_uncertainties(uncertainties_at_ech::UncertaintiesAtEch, injection_name::String, ts::Dates.DateTime)
    injection_uncertainties = get_uncertainties(uncertainties_at_ech, injection_name)
    return get_uncertainties(injection_uncertainties, ts)
end
function get_uncertainties(uncertainties::Uncertainties, ech::Dates.DateTime, injection_name::String, ts::Dates.DateTime)
    uncertainties_at_ech = get_uncertainties(uncertainties, ech)
    return get_uncertainties(uncertainties_at_ech, injection_name, ts)
end

function get_uncertainties(injection_uncertainties::InjectionUncertainties, ts::Dates.DateTime, scenario::String)
    return get_uncertainties(injection_uncertainties, ts)[scenario]
end
function get_uncertainties(uncertainties_at_ech::UncertaintiesAtEch, injection_name::String, ts::Dates.DateTime, scenario::String)
    injection_uncertainties = get_uncertainties(uncertainties_at_ech, injection_name)
    return get_uncertainties(injection_uncertainties, ts, scenario)
end
function get_uncertainties(uncertainties::Uncertainties, ech::Dates.DateTime, injection_name::String, ts::Dates.DateTime, scenario::String)
    uncertainties_at_ech = get_uncertainties(uncertainties, ech)
    return get_uncertainties(uncertainties_at_ech, injection_name, ts, scenario)
end


function sum_uncertainties(uncertainties_at_ech::UncertaintiesAtEch,
    ids::Vector{String},
    ts::Dates.DateTime, s::String)
    sum = 0.
    for id in ids
        sum += get_uncertainties(uncertainties_at_ech, id, ts)[s]
    end
    return sum
end

function compute_prod(uncertainties_at_ech::UncertaintiesAtEch,
    limitables_gens, ts::Dates.DateTime, s::String)
    generators_ids = map(Networks.get_id, limitables_gens)
    return sum_uncertainties(uncertainties_at_ech, generators_ids, ts, s)
end
function compute_prod(uncertainties_at_ech::UncertaintiesAtEch,
    network::Networks.Network, ts::Dates.DateTime, s::String)
    limitables_list = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    return compute_prod(uncertainties_at_ech, limitables_list, ts, s)
end
function compute_prod(uncertainties::Uncertainties,
    network::Networks.Network, ech::Dates.DateTime, ts::Dates.DateTime, s::String)
    uncertainties_at_ech = get_uncertainties(uncertainties, ech)
    return compute_prod(uncertainties_at_ech, network, ts, s)
end

function compute_load(uncertainties_at_ech::UncertaintiesAtEch,
    buses_list, ts::Dates.DateTime, s::String)
    buses_ids = map(Networks.get_id, buses_list)
    return sum_uncertainties(uncertainties_at_ech, buses_ids, ts, s)
end
function compute_load(uncertainties_at_ech::UncertaintiesAtEch,
    network::Networks.Network, ts::Dates.DateTime, s::String)
    return compute_load(uncertainties_at_ech, Networks.get_buses(network), ts, s)
end
function compute_load(uncertainties::Uncertainties,
    network::Networks.Network, ech::Dates.DateTime, ts::Dates.DateTime, s::String)
    uncertainties_at_ech = get_uncertainties(uncertainties, ech)
    return compute_load(uncertainties_at_ech, network, ts, s)
end
