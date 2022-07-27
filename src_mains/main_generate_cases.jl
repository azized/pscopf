"""
    main_ptdf
A main file to generate a PTDF file

Parameters:
    input_path : path to input data directory describing a grid
                (not the pscopf_ files but branches.txt and buses.txt)
"""

using Random
using Dates
using Printf
using StatsBase
using Distributions
using DataStructures
using JuMP
using TimerOutputs

using PSCOPF

root_path = dirname(@__DIR__)
push!(LOAD_PATH, root_path);
cd(root_path)
include(joinpath(root_path, "src", "PTDF.jl"));


###############################################
# Definitions & utils
###############################################

struct PilotableTemplate
    name::String
    p_min::Float64
    p_max::Float64
    start_cost::Float64
    prop_cost::Float64
    dmo::Second
    dp::Second
end


function add_template_to_bus!(network::PSCOPF.Network, bus::PSCOPF.Bus,
                            template::PilotableTemplate)
    unit_name = @sprintf("%s_%s", PSCOPF.get_id(bus), template.name)
    PSCOPF.Networks.add_new_generator_to_bus!(network, PSCOPF.get_id(bus), unit_name, PSCOPF.Networks.PILOTABLE,
                                                template.p_min, template.p_max, #pmin, pmax
                                                template.start_cost, template.prop_cost, #start_cost, prop_cost
                                                template.dmo, template.dp) #dmo, dp
end


"""
probs : list st probs[i] is the probability that a given bus holds exactly i generators
        The probability that a bus holds no generators is 1-sum(probs)
"""
function add_pilotable_generators!(network::PSCOPF.Network, 
                                nb_generators_probabilities::Vector{Float64},
                                # template_probabilities::Vector{Float64},
                                pilotables_templates::Vector{PilotableTemplate})
    @assert (length(nb_generators_probabilities) == length(pilotables_templates))
    @assert all(p>-1e-09 for p in nb_generators_probabilities)
    @assert ( 1 >= sum(nb_generators_probabilities) > 0)

    no_generators_prob = 1 - sum(nb_generators_probabilities)
    nb_generators_probs = ProbabilityWeights(vcat(no_generators_prob, nb_generators_probabilities))
    nb_generators = length(pilotables_templates)

    for bus in PSCOPF.get_buses(network)
        n_generators_on_bus = sample(0:nb_generators, nb_generators_probs)
        templates_to_add = sample(pilotables_templates, n_generators_on_bus, replace=false)
        for pilotable_template in templates_to_add
            add_template_to_bus!(network, bus, pilotable_template)
        end
    end

end


function generate_initial_network(ptdf_network::PTDF.Network,
                                ptdf_folder,
                                default_limit,
                                nb_generators_probabilities, pilotables_templates,
                                )::PSCOPF.Networks.Network
    network = PSCOPF.Networks.Network("generated_network")

    #Buses
    #######
    buses_ids::Set{String} = Set{String}(bus.name for bus in values(ptdf_network.buses))
    PSCOPF.add_new_buses!(network, buses_ids);

    #Branches
    ##########
    branches_ids = Set{String}(branch.name for  branch in values(ptdf_network.branches))
    n_1_cases = branches_ids
    for branch_id in branches_ids
        PSCOPF.add_new_branch!(network, branch_id, default_limit);
        for network_case in n_1_cases
            PSCOPF.add_new_limit!(network, branch_id, network_case, default_limit)
        end
    end

    #PTDF
    ######
    PSCOPF.PSCOPFio.read_ptdf!(network, ptdf_folder)

    #Limitables
    ############
    # None for now
    @warn "No Limitables in the network"

    #Pilotables
    ############
    add_pilotable_generators!(network, nb_generators_probabilities, pilotables_templates)

    return network
end

function generate_init_state(initial_network, output_folder)
    gen_init = SortedDict{String, PSCOPF.GeneratorState}(
        PSCOPF.get_id(gen) => PSCOPF.ON
            for gen in PSCOPF.get_generators(initial_network)
                if PSCOPF.needs_commitment(gen)
    )

    return gen_init
end

"""
randomly distribute a ratio (conso_ratio) of the network's capacity (sum of p_max) on buses
"""
function generate_base_consumptions(network::PSCOPF.Network,
                                    ech::DateTime, ts::DateTime, base_s::String,
                                    conso_ratio::Float64)::PSCOPF.Uncertainties
    @assert ( 0 < conso_ratio <= 1. )

    uncertainties = PSCOPF.Uncertainties()

    network_capacity = sum( PSCOPF.get_p_max(generator)
                            for generator in PSCOPF.get_generators_of_type(network, PSCOPF.PILOTABLE))
    total_consumptio = network_capacity * conso_ratio

    nb_buses = PSCOPF.get_nb_buses(network)
    buses_ids = PSCOPF.get_id.(PSCOPF.get_buses(network))
    bus_conso_ratios = rand(Dirichlet(nb_buses, 1.0))
    distribution = Dict(zip(buses_ids, bus_conso_ratios))

    for (bus_id, bus_conso_ratio) in distribution
        conso_l = bus_conso_ratio * total_consumptio
        PSCOPF.add_uncertainty!(uncertainties, ech, bus_id, ts, base_s, conso_l)
    end

    return uncertainties
end

# uncertainties[ech][nodal_injection_name][ts][scenario_name]
function multiply(uncertainties::PSCOPF.Uncertainties, coeff::Number)
    result_uncertainties = PSCOPF.Uncertainties()
    for (ech, uncertainties_at_ech) in uncertainties
        result_uncertainties_at_ech = get!(result_uncertainties, ech, PSCOPF.UncertaintiesAtEch())
        for (inj_name, by_ts_uncerts) in uncertainties_at_ech
            result_by_ts_uncerts = get!(result_uncertainties_at_ech, inj_name, PSCOPF.InjectionUncertainties())
            for (ts, by_scenario_uncerts) in by_ts_uncerts
                result_by_scenario_uncerts = get!(result_by_ts_uncerts, ts, SortedDict{String, Float64}())
                for (s, val) in by_scenario_uncerts
                    result_by_scenario_uncerts[s] = coeff*val
                end
            end
        end
    end
    return result_uncertainties
end


function compute_free_flows(network::PSCOPF.Network, ech, ts, gen_init, uncertainties, output_folder)
    custom_mode = PSCOPF.ManagementMode("custom_mode", Dates.Minute(60))
    tso = PSCOPF.TSOOutFO(PSCOPF.TSOConfigs(CONSIDER_N_1_CSTRS=true))
    initial_context = PSCOPF.PSCOPFContext(network, [ts], custom_mode,
                                    gen_init,
                                    uncertainties, nothing,
                                    output_folder)
    result, _ = PSCOPF.run_step!(initial_context, tso, ech, nothing)

    free_flows = SortedDict{Tuple{String,String},Float64}()
    for ((branch_id,_,_,ptdf_case), flow_expr) in result.flows
        free_flows[branch_id, ptdf_case] = value(flow_expr)
    end

    return free_flows
end


function update_network_limits!(network::PSCOPF.Network, flows::SortedDict{Tuple{String,String},Float64}, ratio)
    for ( (branch_id, network_case), flow_l) in flows
        limit_l = ceil(abs(flow_l * ratio))
        PSCOPF.add_new_limit!(network, branch_id, network_case, limit_l)
    end
end


function generate_uncertainties(network, base_values::PSCOPF.UncertaintiesAtEch, ts, base_s, ECH, nb_scenarios, prediction_error)
    TS = PSCOPF.create_target_timepoints(ts)
    uncertainties_distributions = SortedDict{String, PSCOPF.UncertaintyErrorNDistribution}()
    for bus in PSCOPF.get_buses(network)
        bus_id = PSCOPF.get_id(bus)
        uncertainties_distributions[bus_id] = PSCOPF.UncertaintyErrorNDistribution(
                                                            bus_id,
                                                            0., 5000,
                                                            PSCOPF.get_uncertainties(base_values, bus_id, ts, base_s),
                                                            prediction_error,
                                                            )
    end

    return PSCOPF.generate_uncertainties(network,
                                        TS, ECH,
                                        uncertainties_distributions,
                                        nb_scenarios)
end


function press_to_continue()
    disable_timer!(PSCOPF.TIMER_TRACKS)
    println("\n"^3)
    println("press enter to continue")
    readline()
    println("\n"^3)
    enable_timer!(PSCOPF.TIMER_TRACKS)
end

"""
1- creates a symlink to the ptdf file if it exists otherwise it computes it (NOTE : existence verification is based on file name not on file content)
2- Generate initial network : randomly generated pilotable units and default limits on branches
3- Generate a base consumption for buses (a single scenario, a single ts)
4- Launch a TSO model to compute base flows using the base consumptions
5- Update Network branch limits using the computed flows and write the partial instance
6- Generate consumption uncertainties randomly based on the base consumptions
"""
function main_instance_generate(input_path,
            ref_bus_num, distributed,
            default_limit, nb_generators_probabilities, pilotables_templates,
            limits_conso_to_unit_capa_ratio,
            limit_to_free_flow_ratio,
            ECH, ts, nb_scenarios, prediction_error,
            output_folder)

    ech = ECH[1]

    ###############################################
    # COMPUTE PTDF for N and N-1
    ###############################################
    @info "#### COMPUTE PTDF N and N-1 ####"
    time_ptdf = @elapsed begin
        ptdf_network = PTDF.read_network(input_path)
        out_ptdf_file = joinpath(output_folder, PTDF.PTDF_FILENAME)
        in_possible_ptdf_file = joinpath(input_path, PTDF.PTDF_FILENAME)
        if isfile(out_ptdf_file)
            @info @sprintf("a PTDF file was already found in %s. PTDF was not recomputed!", out_ptdf_file)
        elseif isfile(in_possible_ptdf_file)
            mkpath(output_folder)
            symlink(abspath(in_possible_ptdf_file), out_ptdf_file)
            @info @sprintf("a PTDF file was already found in %s. Created a symlink to the existing file!", in_possible_ptdf_file)
        else
            PTDF.compute_and_write_n_non_bridges(ptdf_network, ref_bus_num, distributed, 1e-6, input_path, output_folder)
        end
    end


    ###############################################
    # Generate Initial Network
    ###############################################
    @info "#### Generate Initial Network ####"
    initial_network = generate_initial_network(ptdf_network,
                                            output_folder,
                                            default_limit,
                                            nb_generators_probabilities, pilotables_templates)
    gen_init = generate_init_state(initial_network, output_folder)
    PSCOPF.PSCOPFio.write(joinpath(output_folder, "initial_network"), initial_network, ignore_ptdf=true)
    symlink(abspath(out_ptdf_file), joinpath(output_folder, "initial_network", PTDF.PTDF_FILENAME))
    PSCOPF.PSCOPFio.write(joinpath(output_folder, "initial_network"), gen_init)


    ###############################################
    # Generate Base Uncertainties for limits
    ###############################################
    @info "#### Generate Base Uncertainties for limits ####"
    base_consumptions = generate_base_consumptions(initial_network, ech, ts, "BASE_S", limits_conso_to_unit_capa_ratio)
    PSCOPF.PSCOPFio.write(joinpath(output_folder, "initial_network"), base_consumptions)


    ###############################################
    # Compute Base Flows
    ###############################################
    @info "#### Compute Base Flows ####"
    time_free_flows = @elapsed free_flows = compute_free_flows(initial_network, ech, ts, gen_init, base_consumptions,
                                                            joinpath(output_folder, "init_limits"))


    instance_path = joinpath(output_folder, "instance")
    ###############################################
    # Generate Network : update branch limits
    ###############################################
    @info "#### Update Limits and Write Network instance ####"
    update_network_limits!(initial_network, free_flows, limit_to_free_flow_ratio)
    generated_network = initial_network
    PSCOPF.PSCOPFio.write(instance_path, generated_network, ignore_ptdf=true)
    symlink(abspath(out_ptdf_file), joinpath(instance_path, PTDF.PTDF_FILENAME))
    PSCOPF.PSCOPFio.write(joinpath(output_folder, "instance"), gen_init)


    ###############################################
    # Generate Uncertainties
    ###############################################
    @info "#### Generate Instance Uncertainties ####"
    base_uncertainties = multiply(base_consumptions, 0.7)
    uncertainties = generate_uncertainties(generated_network, base_uncertainties[ech], ts, "BASE_S", ECH, nb_scenarios, prediction_error)
    PSCOPF.PSCOPFio.write(instance_path, uncertainties)

    return (generated_network, gen_init, uncertainties), (time_ptdf, time_free_flows)
end

###############################################
# INPUT & PARAMS
###############################################
Random.seed!(0)

matpower_case = "case5"
input_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "data_matpower", matpower_case) )
output_folder = joinpath(@__DIR__, "..", "data", matpower_case)

# PTDF
#######
ref_bus_num = 1
distributed = true

# Initial Network
##################
default_limit = 1e5
pilotables_templates = [
    PilotableTemplate("_15m",   0., 500.,    0., 30.,  Second(Minute(15)), Second(Minute(15)))
    PilotableTemplate("_.5h",  50., 200.,  300., 27.,  Second(Minute(30)), Second(Minute(15)))
    PilotableTemplate("_1h",   50., 200.,  250., 25.,  Second(Hour(1)),    Second(Minute(15)))
    PilotableTemplate("_2h",   50., 300.,  100., 20.,  Second(Hour(2)),    Second(Minute(15)))
    PilotableTemplate("_4h",   50., 600.,  150., 15.,  Second(Hour(4)),    Second(Minute(15)))
]
nb_generators_probabilities = [.25, .2, .25, .1, .05] #no_generator_proba : 0.15
@assert (length(nb_generators_probabilities) == length(pilotables_templates))


# Base Uncertainties
#####################
limits_conso_to_unit_capa_ratio = 1. #consumption volumne is equal to maximum generators capacity (all generators produce at max capacity : optim is useless !)
# limits_conso_to_unit_capa_ratio = 0.7 #consumption used for branch dimensioning will represent ?% of the units' max capacities (distributed randomly)


# Limits
#########
limit_to_free_flow_ratio = 1.

# Uncertainties
################
# base_uncertainties = base_uncertainties_coeff * limits_conso_to_unit_capa_ratio * maxNetworkCapacity
# (keeping same distribution of the base consumption used to compute limits)
base_uncertainties_coeff = 0.7
ts1 = DateTime("2015-01-01T11:00:00")
ECH = [ts1-Hour(4), ts1-Hour(2), ts1-Hour(1), ts1-Minute(30), ts1-Minute(15), ts1]
nb_scenarios = 3
prediction_error = 0.01


#################################################################################################################
# Launch
#################################################################################################################
logfile = PSCOPF.get_config("TEMP_GLOBAL_LOGFILE")
open(logfile, "a") do file_l
    write(file_l, "-"^120 * "\n")
    write(file_l, @sprintf("generation of usecase : %s\n", output_folder))
    write(file_l, @sprintf("n-1? : %s\n", PSCOPF.get_config("CONSIDER_N_1")))
end


time_generation = @elapsed (generated_network, gen_init, uncertainties), (time_ptdf, time_free_flows) =
                            main_instance_generate(input_path,
                                ref_bus_num, distributed,
                                default_limit, nb_generators_probabilities, pilotables_templates,
                                limits_conso_to_unit_capa_ratio,
                                limit_to_free_flow_ratio,
                                ECH, ts1, nb_scenarios, prediction_error,
                                output_folder)

println("Computing all ptdfs took:", time_ptdf)
println("Computing free flows took :", time_free_flows)
println("Generating the whole network instance took :", time_generation)

println(PSCOPF.TIMER_TRACKS)

nb_rso_constraint = PSCOPF.nb_rso_constraint(generated_network, nb_scenarios, 4)
open(logfile, "a") do file_l
    write(file_l, @sprintf("nb rso constraints : %d\n", nb_rso_constraint))
    write(file_l, @sprintf("Computing all ptdfs took: %s\n", time_ptdf))
    write(file_l, @sprintf("Computing free flows took : %s\n", time_free_flows))
    write(file_l, @sprintf("Generating the whole network instance took : %s\n", time_generation))
    write(file_l, @sprintf("TIMES:\n%s\n", PSCOPF.TIMER_TRACKS))
end
