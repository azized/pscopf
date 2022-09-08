"""
    main_illustration
A main file to launch PSCOPF illustrating the cases in data/illustrations
"""

using Dates
using DataStructures

root_path = dirname(@__DIR__)
push!(LOAD_PATH, root_path);
cd(root_path)
include(joinpath(root_path, "src", "PSCOPF.jl"));

PSCOPF.set_config!("LP_FILES", true)


function read_context(instance_path, outname)
    output_path = joinpath(instance_path, outname)
    PSCOPF.rm_non_prefixed(output_path, "pscopf_")

    TS = [Dates.DateTime("2015-01-01T11:00:00")]
    network = PSCOPF.Data.pscopfdata2network(instance_path)
    uncertainties = PSCOPF.PSCOPFio.read_uncertainties(instance_path)
    generators_init_state = PSCOPF.PSCOPFio.read_initial_state(instance_path)

    exec_context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.ManagementMode("custom", Dates.Minute(1)),
                                        generators_init_state,
                                        uncertainties, nothing,
                                        output_path)
    return exec_context
end

function usecase_1_1(instance_path)
    ECH = [Dates.DateTime("2015-01-01T07:00:00")]
    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOOutFO()],
        ))

    # load Data
    exec_context = read_context(instance_path, "tso")

    PSCOPF.run!(exec_context, sequence)
end
function usecase_1_2(instance_path)
    ECH = [Dates.DateTime("2015-01-01T07:00:00")]
    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOBilevel()],
        ))

    # load Data
    exec_context = read_context(instance_path, "tsobilevel")

    PSCOPF.run!(exec_context, sequence)
end
function usecase_1_3(instance_path)
    ECH = [Dates.DateTime("2015-01-01T07:00:00")]
    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(CONSIDER_DELTAS=false))],
        ))

    # load Data
    exec_context = read_context(instance_path, "tsobilevel_nodeltas")

    PSCOPF.run!(exec_context, sequence)
end




#########################
# EXECUTION
#########################
IN = joinpath(@__DIR__, "..", "data", "illustrations", "usecase1")
usecase_1_3(IN)
