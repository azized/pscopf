using Distributions
using Dates

using .Networks

struct UncertaintiesGenerator <: AbstractDataGenerator
    network::Networks.Network
    target_timepoints::Vector{Dates.DateTime}
    horizon_timepoints::Vector{Dates.DateTime}
    #potentially split into load_uncertainties_distribution and unit_uncertainties_distribution :
    uncertainties_distribution #ditribution info for each limitable unit and load
    nb_scenarios::Int
end

function scenario_name(index_scenario::Int)
    return "S"*string(index_scenario)
end

function generate_uncertainties(network::Networks.Network,
                                target_timepoints::Vector{Dates.DateTime},
                                horizon_timepoints::Vector{Dates.DateTime},
                                uncertainties_distribution, nb_scenarios::Int)
    generator = UncertaintiesGenerator(network, target_timepoints, horizon_timepoints, uncertainties_distribution, nb_scenarios)
    return launch(generator)
end

function launch(uncertainties_generator::UncertaintiesGenerator)
    uncertainties = Uncertainties()
    if !check_uncertainties_distribution(uncertainties_generator.network, uncertainties_generator.uncertainties_distribution)
        error("Invalid uncertainties parameters!")
    end

    scenarios = scenario_name.(1:uncertainties_generator.nb_scenarios)

    for ech in uncertainties_generator.horizon_timepoints
        for bus_or_gen in union(Networks.get_buses(uncertainties_generator.network),
                                Networks.get_generators_of_type(uncertainties_generator.network, Networks.LIMITABLE))
            nodal_injection_name = Networks.get_id(bus_or_gen)
            uncertain_distro = uncertainties_generator.uncertainties_distribution[nodal_injection_name]
            for ts in uncertainties_generator.target_timepoints
                values = generate_values(uncertain_distro, ech, ts, uncertainties_generator.nb_scenarios)
                for (scenario_name,val) in zip(scenarios,values)
                    add_uncertainty!(uncertainties, ech, nodal_injection_name, ts, scenario_name, val)
                end
            end
        end
    end

    return uncertainties
end

function generate_values(uncertain_distro::UncertaintyNDistribution,
                        ech::Dates.DateTime, ts::Dates.DateTime,
                        nb_scenarios::Int64)

    delta_time_l = max(0,
                    Dates.value(floor(ts - ech, Dates.Second)) / 3600)
    factor_l = delta_time_l * uncertain_distro.time_factor
    adjusted_sigma_l = factor_l * uncertain_distro.sigma
    random_values_l = rand(Distributions.Normal(uncertain_distro.mu, adjusted_sigma_l), nb_scenarios)
    values_l = max.(uncertain_distro.min_value, random_values_l)
    values_l = min.(uncertain_distro.max_value, values_l)

    return values_l
end

function generate_values(uncertain_distro::UncertaintyErrorNDistribution,
                        ech::Dates.DateTime, ts::Dates.DateTime,
                        nb_scenarios::Int64)

    delta_time_l = max(0,
                    Dates.value(floor(ts - ech, Dates.Second)) / 3600)
    adjusted_sigma_l = delta_time_l * uncertain_distro.error_sigma
    rand_deviations_l = rand(Distributions.Normal(0., adjusted_sigma_l), nb_scenarios)
    random_values_l = uncertain_distro.mu * (1 .+ rand_deviations_l)
    values_l = max.(uncertain_distro.min_value, random_values_l)
    values_l = min.(uncertain_distro.max_value, values_l)

    return values_l
end

############################################################
#         Checkers
############################################################

function check_extra_rows(network, uncertainties_distribution)
    #warn if extra input
    for (id, _) in uncertainties_distribution
        gen_or_bus = Networks.get_generator_or_bus(network, id)
        if ( ismissing(gen_or_bus) ||
            ( isa(gen_or_bus, Networks.Generator) && (Networks.get_type(gen_or_bus) != Networks.LIMITABLE) )
            )
            msg = @sprintf("Uncertainty distribution input for id %s may not be needed.", id)
            @warn(msg)
        end
    end

    return true
end

function check_buses_listed(network, uncertainties_distribution)
    checks = true

    #all buses are listed
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        if !haskey(uncertainties_distribution, bus_id)
            msg = @sprintf("Missing uncertainty distribution for bus %s.", bus_id)
            @error(msg)
            checks = false
        end
    end

    return checks
end


function check_limitables_listed(network, uncertainties_distribution)
    checks = true

    #all limitables are listed
    for gen in Networks.get_generators_of_type(network, Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        if !haskey(uncertainties_distribution, gen_id)
            msg = @sprintf("Missing uncertainty distribution for limitable %s.", gen_id)
            @error(msg)
            checks = false
        end
    end

    return checks
end

function check_min_max(network, uncertainties_distribution)
    checks = true

    #max_value of limitables are respected
    for (id, uncertain_distro) in uncertainties_distribution
        @assert(id == uncertain_distro.id)
        gen_or_bus = Networks.get_generator_or_bus(network, id)
        if !ismissing(gen_or_bus) && isa(gen_or_bus, Networks.Generator)
            if Networks.get_p_max(gen_or_bus) < uncertain_distro.max_value
                msg = @sprintf("Invalid Uncertainty distribution max value (%f) for id %s. \
                                Value must be less or equal to %f (the generator's Pmax).",
                                uncertain_distro.max_value, id, Networks.get_p_max(gen_or_bus))
                @error(msg)
                checks = false
            end
        end
    end

    #  max_value >= min_value >= 0
    for (id, uncertain_distro) in uncertainties_distribution
        if uncertain_distro.max_value < uncertain_distro.min_value
            msg = @sprintf("Invalid Uncertainty distribution for id %s. \
                                max_value (%f) must be greater or equal to min_value (%f).",
                                id, uncertain_distro.max_value, uncertain_distro.min_value)
                @error(msg)
                checks = false
        end
        if uncertain_distro.min_value < 0
            msg = @sprintf("Invalid Uncertainty distribution for id %s. \
                                pmin (%f) must be non-negative.",
                                id, uncertain_distro.min_value)
                @error(msg)
                checks = false
        end
    end

    return checks
end

function check_uncertainties_distribution(network, uncertainties_distribution)
    return ( check_extra_rows(network, uncertainties_distribution)
            & check_buses_listed(network, uncertainties_distribution)
            & check_limitables_listed(network, uncertainties_distribution))
end

function check_uncertainties_distribution(network, uncertainties_distribution::UncertaintiesDistribution)
    return ( check_extra_rows(network, uncertainties_distribution)
            & check_buses_listed(network, uncertainties_distribution)
            & check_limitables_listed(network, uncertainties_distribution)
            & check_min_max(network, uncertainties_distribution) )
end
