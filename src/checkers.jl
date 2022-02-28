using Dates
using DataStructures

using ..Networks

###########################################################
###         Network Checkers
###########################################################

function check(network::Network)
    checks = true

    for generator in Networks.get_generators(network)
        checks &= check(generator)
    end

    for branch in Networks.get_branches(network)
        checks &= check(branch)
    end

    checks &= check_ptdf(network.ptdf, network)

    return checks
end

####    PTDF
##################################

function check_ptdf_branch_entries(ptdf::SortedDict{String,SortedDict{String, Float64}},
                                    network::Network)
    checks = true

    network_branch_ids = Set{String}(map(Networks.get_id, Networks.get_branches(network)))
    ptdf_branch_ids = Set{String}(keys(ptdf))

    missing_entries = setdiff(network_branch_ids, ptdf_branch_ids)
    if !isempty(missing_entries)
        checks = false
        msg = @sprintf("Missing branch entries `%s` in ptdf", missing_entries)
        @error(msg)
    end

    extra_entries = setdiff(ptdf_branch_ids, network_branch_ids)
    if !isempty(extra_entries)
        checks = false
        msg = @sprintf("Extra branch entries `%s` in ptdf", extra_entries)
        @error(msg)
    end

    return checks
end
function check_ptdf_bus_entries(ptdf::SortedDict{String,SortedDict{String, Float64}},
                                network::Network)
    checks = true
    network_bus_ids = Set{String}(map(Networks.get_id, Networks.get_buses(network)))
    for (branch_id, branch_ptdf) in ptdf
        ptdf_bus_ids = Set{String}(keys(branch_ptdf))

        missing_entries = setdiff(network_bus_ids, ptdf_bus_ids)
        if !isempty(missing_entries)
            checks = false
            msg = @sprintf("Missing bus entries `%s` for branch `%s` in ptdf", missing_entries, branch_id)
            @error(msg)
        end

        extra_entries = setdiff(ptdf_bus_ids, network_bus_ids)
        if !isempty(extra_entries)
            checks = false
            msg = @sprintf("Extra bus entries `%s` for branch `%s` in ptdf", extra_entries, branch_id)
            @error(msg)
        end

    end
    return checks
end
function check_ptdf(ptdf::SortedDict{String,SortedDict{String, Float64}},
                    network::Network)
    return ( check_ptdf_branch_entries(ptdf, network)
            && check_ptdf_bus_entries(ptdf, network) )
end

####    Generator
##################################

function check(generator::Generator)
    checks = true

    msg_prefix = @sprintf("Generator %s : ", Networks.get_id(generator))

    if Networks.get_p_min(generator) < 0
        checks = false
        msg = @sprintf("Invalid input %f : p_min must be non-negative.", Networks.get_p_min(generator))
        @error(msg_prefix * msg)
    end

    if Networks.get_p_max(generator) < Networks.get_p_min(generator)
        checks = false
        msg = @sprintf("Invalid input %f : p_max must be greater than or equal to p_min (i.e. %f).", Networks.get_p_max(generator), Networks.get_p_min(generator))
        @error(msg_prefix * msg)
    end

    if Networks.get_start_cost(generator) < 0
        checks = false
        msg = @sprintf("Invalid input %f : start_cost must be non-negative.", Networks.get_start_cost(generator))
        @error(msg_prefix * msg)
    end

    if Networks.get_prop_cost(generator) < 0
        checks = false
        msg = @sprintf("Invalid input %f : prop_cost must be non-negative.", Networks.get_prop_cost(generator))
        @error(msg_prefix * msg)
    end

    if Networks.get_dp(generator) < Dates.Second(0)
        checks = false
        msg = @sprintf("Invalid input %s : dp must be non-negative.", Networks.get_dp(generator))
        @error(msg_prefix * msg)
    end

    if Networks.get_dmo(generator) < Networks.get_dp(generator)
        checks = false
        msg = @sprintf("Invalid input %s : dmo must be greater than or equal to dp (i.e. %s).",
                        Networks.get_dmo(generator), Networks.get_dp(generator))
        @error(msg_prefix * msg)
    end

    #if pmin=0, we can suppose that the generator is always
    # 1- the start cost is supposed to be paid in the far past => 0
    # 2- no need for DMO since the decision is always ON => DMO=DP (actually, infinite DMO)
    if Networks.get_p_min(generator) < 1e-09
        if Networks.get_start_cost(generator) > 0
            checks = false
            msg = @sprintf("Invalid input %s : generators with null p_min must have 0 start cost.", Networks.get_start_cost(generator))
            @error(msg_prefix * msg)
        end

        if Networks.get_dmo(generator) != Networks.get_dp(generator)
            checks = false
            msg = @sprintf("Invalid input %s : generators with null p_min must have DMO=DP (%s).",
                        Networks.get_dmo(generator), Networks.get_dp(generator))
            @error(msg_prefix * msg)
        end
    end

    if Networks.get_type(generator) == Networks.LIMITABLE
        if Networks.get_p_min(generator) > 0
            checks = false
            msg = @sprintf("Invalid input %f : Limitable units must have a minimum production capacity of 0.", Networks.get_p_min(generator))
            @error(msg_prefix * msg)
        end
    end

    return checks
end

####    Branch
##################################

function check(branch::Branch)
    if Networks.get_limit(branch) < 0
        msg = @sprintf("Branch %s: Invalid input %f : limit must be non-negative.",
                    Networks.get_id(branch), Networks.get_limit(branch))
        @error(msg)
        return false
    end
    return true
end


###########################################################
###         Uncertainties
###########################################################

function check_uncertainties(uncertainties::Uncertainties, network)
    return ( check_uncertainties_same_scenarios(uncertainties)
            && check_uncertainties_same_timesteps(uncertainties)
            && check_uncertainties_values(uncertainties, network)
            && check_uncertainties_underprod(uncertainties, network)
            && check_uncertainties_limitables(uncertainties, network)
            && check_uncertainties_buses(uncertainties, network)
    )
end

function check_uncertainties_same_scenarios(uncertainties::Uncertainties)
    reference_scenarios = get_scenarios(uncertainties)
    for (ech,_) in uncertainties
        for (injection_name,_) in get_uncertainties(uncertainties, ech)
            for (ts, by_scenario_injections) in get_uncertainties(uncertainties, ech, injection_name)
                scenarios = collect(keys(by_scenario_injections))
                if scenarios != reference_scenarios
                    msg = @sprintf("Different scenarios list at ech=%s for name=%s and ts=%s. (Reference scenarios : %s)",
                                    ech, injection_name, ts, reference_scenarios)
                    @error(msg)
                    return false
                end
            end
        end
    end
    return true
end

function check_uncertainties_same_timesteps(uncertainties::Uncertainties)
    reference_TS = get_target_timepoints(uncertainties)
    for (ech,_) in uncertainties
        for (injection_name, by_ts_uncertainties) in get_uncertainties(uncertainties, ech)
            TS = collect(keys(by_ts_uncertainties))
            if TS != reference_TS
                msg = @sprintf("Different timesteps list at ech=%s for name=%s. (reference TS : %s)",
                                ech, injection_name, reference_TS)
                @error(msg)
                return false
            end
        end
    end
    return true
end


function check_uncertainties_values(uncertainties::Uncertainties, network::Network)
    checks = true
    for (ech,_) in uncertainties
        for (injection_name,_) in get_uncertainties(uncertainties, ech)
            gen_or_bus = Networks.get_generator_or_bus(network, injection_name)
            if ( ismissing(gen_or_bus)
                || (typeof(gen_or_bus) <: Networks.Generator) && (Networks.get_type(gen_or_bus) != Networks.LIMITABLE) )
                msg = @sprintf("Invalid injection name `%s` : injection names must refer to buses or limitables ids",
                                injection_name)
                @error(msg)
                checks = false
                continue
            end

            for (ts, _) in get_uncertainties(uncertainties, ech, injection_name)
                for (scenario, value) in get_uncertainties(uncertainties, ech, injection_name, ts)
                    if value < 0
                        msg = @sprintf("Invalid injection value %f at (ech=%s, name=%s, ts=%s, s=%s) : injections must be non-negative",
                                    value, ech, injection_name, ts, scenario)
                        @error(msg)
                        checks = false
                    elseif ( ( typeof(gen_or_bus) <: Networks.Generator )
                            && value > Networks.get_p_max(gen_or_bus) )
                        msg = @sprintf("Invalid injection value %f at (ech=%s, name=%s, ts=%s, s=%s) : injection is greater than pmax (%f)",
                        value, ech, injection_name, ts, scenario, Networks.get_p_max(gen_or_bus))
                        @error(msg)
                        checks = false
                    end
                end
            end
        end
    end
    return checks
end

"""
    check that for each ech, ts, s:
        prod < load
    To be sure that we can automatically set limitables to the available production levels
    without exceeding demand
    and then compliment missing demand with imposables
"""
function check_uncertainties_underprod(uncertainties::Uncertainties, network)
    checks = true
    target_timepoints = get_target_timepoints(uncertainties)
    scenarios = get_scenarios(uncertainties)

    for (ech,uncertainties_at_ech) in uncertainties
        for ts in target_timepoints
            for scenario in scenarios
                load = compute_load(uncertainties_at_ech, network, ts, scenario)
                prod = compute_prod(uncertainties_at_ech, network, ts, scenario)
                if prod > load
                    msg = @sprintf("Invalid injections at (ech=%s, ts=%s, s=%s) : \
                                    renewables injections (%f) must not be greater than load (%f) !",
                                    ech, ts, scenario, prod, load)
                    @error(msg)
                    checks = false
                end
            end
        end
    end
    return checks
end

"""
    make sure all limitables have listed values
"""
function check_uncertainties_limitables(uncertainties::Uncertainties, network)
    checks = true
    limitables = Networks.get_generators_of_type(network, Networks.LIMITABLE)
    limitables_ids = Set{String}( map(Networks.get_id, limitables) )
    for (ech,uncertainties_at_ech) in uncertainties
        injection_ids = Set{String}( keys(uncertainties_at_ech) )
        for missing_limitable_id in setdiff(limitables_ids, injection_ids)
            msg = @sprintf("Missing value for limitable injection `%s` at ech=%s",
                            missing_limitable_id, ech)
            @error(msg)
            checks = false
        end
    end
    return checks
end

"""
    make sure all buses have listed values
"""
function check_uncertainties_buses(uncertainties::Uncertainties, network)
    checks = true
    buses_ids = Set{String}( map(Networks.get_id, Networks.get_buses(network)) )
    for (ech,uncertainties_at_ech) in uncertainties
        injection_ids = Set{String}( keys(uncertainties_at_ech) )
        for missing_bus_id in setdiff(buses_ids, injection_ids)
            msg = @sprintf("Missing value for load `%s` at ech=%s",
                            missing_bus_id, ech)
            @error(msg)
            checks = false
        end
    end
    return checks
end


###########################################################
###         Context
###########################################################

function check(context)
    return (
        check(get_network(context))
        && check_uncertainties(get_uncertainties(context), get_network(context))
        && check_initial_state(get_generators_initial_state(context), get_network(context))
        && check_target_timepoints(get_target_timepoints(context))
        && check_uncertainties_contain_ts(get_uncertainties(context), get_target_timepoints(context))
        && check_uncertainties_contains_ech(get_uncertainties(context), get_horizon_timepoints(context))
    )
end

####    Initial State
##################################

function check_initial_state(initial_state::SortedDict{String, GeneratorState}, network)
    checks = true
    must_list_generators = filter(gen -> Networks.get_p_min(gen) > 0,
                                collect(Networks.get_generators(network)) )
    must_list_gen_ids = Set(map( gen -> Networks.get_id(gen),
                                must_list_generators ))
    listed_gen_ids = Set(keys(initial_state))

    for missing_gen_id in setdiff(must_list_gen_ids, listed_gen_ids)
        msg = @sprintf("Initial state for generator %s is missing", missing_gen_id)
        @error(msg)
        checks = false
    end

    for extra_gen_id in setdiff(listed_gen_ids, must_list_gen_ids)
        if initial_state[extra_gen_id] == OFF
            msg = @sprintf("Initial state for generator %s must be ON (or can be ommited). \
                            Generators with no Pmin can always be considered on and no start cost is linked to them",
                            extra_gen_id)
            @error(msg)
            checks = false
        end
    end

    return checks
end

####    TS
##################################

function check_target_timepoints(target_timepoints)
    if !issorted(target_timepoints)
        @error("target timepoints vector is not sorted")
        return false
    end

    if length(target_timepoints) > length(unique(target_timepoints))
        @error("target timepoints vector contains duplicate values")
        return false
    end

    return true
end


####    Uncertainties-related
##################################

function check_uncertainties_contain_ts(uncertainties::Uncertainties, target_timepoints)
    checks = true
    uncertainties_target_timepoints = get_target_timepoints(uncertainties)
    if !issubset(target_timepoints, uncertainties_target_timepoints)
        for ts in setdiff(target_timepoints, uncertainties_target_timepoints)
            msg = @sprintf("Target timepoint %s is missing in uncertainties", ts)
            @error(msg)
            checks = false
        end
    end
    return checks
end

function check_uncertainties_contains_ech(uncertainties::Uncertainties, horizon_timepoints)
    checks = true
    uncertainties_horizon_timepoints = get_horizon_timepoints(uncertainties)
    if !issubset(horizon_timepoints, uncertainties_horizon_timepoints)
        for ech in setdiff(horizon_timepoints, uncertainties_horizon_timepoints)
            msg = @sprintf("Horizon timepoint %s is missing in uncertainties", ech)
            @error(msg)
            checks = false
        end
    end
    return checks
end
