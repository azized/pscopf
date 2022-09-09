"""
    main_demo
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

    network = PSCOPF.Data.pscopfdata2network(instance_path)
    uncertainties = PSCOPF.PSCOPFio.read_uncertainties(instance_path)
    generators_init_state = PSCOPF.PSCOPFio.read_initial_state(instance_path)

    ts1 = Dates.DateTime("2015-01-01T11:00:00") #11h
    TS = PSCOPF.create_target_timepoints(ts1) #T: 11h, 11h15, 11h30, 11h45
    ECH = PSCOPF.generate_ech(network, TS, PSCOPF.PSCOPF_MODE_1) #ech: -4h, -1h, -30mins, -15mins, 0h

    exec_context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.ManagementMode("custom", Dates.Minute(1)),
                                        generators_init_state,
                                        uncertainties, nothing,
                                        output_path)
    return exec_context, ECH
end

function usecase_EOD(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_eod")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_RSO(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_rso")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.TSOOutFO()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_seq1(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_seq1")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOOutFO()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_seq15(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_seq15")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.TSOOutFO(), PSCOPF.EnergyMarket()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_seq2(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_seq2")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOBilevel(), PSCOPF.BalanceMarket()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_link(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_lttd")

    sequence = PSCOPF.Sequence(SortedDict(
            ECH[1] => [PSCOPF.EnergyMarket(),
                        PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL=true)),
                        PSCOPF.BalanceMarket()],
            ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_mode1(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_mode1")

    sequence = PSCOPF.Sequence(SortedDict(
        ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOOutFO()],
        ECH[2] => [PSCOPF.EnergyMarketAtFO(), PSCOPF.TSOOutFO()],
        ECH[3] => [PSCOPF.TSOOutFO()],
        ECH[4] => [PSCOPF.TSOOutFO()],
        ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end

function usecase_mode2(instance_path)
    # load Data
    exec_context, ECH = read_context(instance_path, "usecase_mode2")

    sequence = PSCOPF.Sequence(SortedDict(
        ECH[1] => [PSCOPF.EnergyMarket(), PSCOPF.TSOOutFO()],
        ECH[2] => [PSCOPF.EnergyMarketAtFO(), PSCOPF.TSOOutFO()],
        ECH[3] => [PSCOPF.EnergyMarket(), PSCOPF.TSOBilevel(), PSCOPF.BalanceMarket()],
        ECH[4] => [PSCOPF.EnergyMarket(), PSCOPF.TSOBilevel(), PSCOPF.BalanceMarket()],
        ECH[end] => [PSCOPF.Assessment()],
        ))

    PSCOPF.run!(exec_context, sequence)
end



#########################
# EXECUTION
#########################
IN = joinpath(@__DIR__, "..", "data", "usecase_demo", "instance")
# usecase_EOD(IN)
# usecase_RSO(IN)
# usecase_seq1(IN)
# usecase_seq15(IN)
# usecase_seq2(IN)
# usecase_link(IN)
# usecase_mode1(IN)
# usecase_mode2(IN)
