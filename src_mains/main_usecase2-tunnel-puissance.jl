"""
    main_mode1
A main file to launch PSCOPF in MODE 1.

Parameters:
    instance_path : path to input data directory
    output_path : path to output directory
"""

using Dates

root_path = dirname(@__DIR__)
push!(LOAD_PATH, root_path);
cd(root_path)
include(joinpath(root_path, "src", "PSCOPF.jl"));




#########################
# INPUT & PARAMS
#########################

# instance_path is the path to the input data directory containing :
# pscopf_units : the description of the units
# pscopf_init : the description of the initial state of units (just before the first interest timepoit T)
# pscopf_limits : the flow capacity of each branch
# pscopf_ptdf : the ptdf coefficients per (branch, bus_id)
# pscopf_uncertainties : the nodal injections (for each bus and each limitable)
instance_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "usecases-euro-simple", "usecase2-tunnel-puissance", "data") )

# output_path is the path where output files will be write_commitment_schedule
#NOTE: all files in output_path, except those starting with pscopf_, will be deleted
output_path = length(ARGS) > 1 ? ARGS[2] : joinpath(instance_path, "..", "output")




#########################
# EXECUTION
#########################

# load Data
network = PSCOPF.Data.pscopfdata2network(instance_path)
uncertainties = PSCOPF.PSCOPFio.read_uncertainties(instance_path)
generators_init_state = PSCOPF.PSCOPFio.read_initial_state(instance_path)

# Launch mode 1
mode = PSCOPF.PSCOPF_MODE_1

ts1 = Dates.DateTime("2015-01-01T11:00:00") #11h
TS = PSCOPF.create_target_timepoints(ts1) #T: 11h, 11h15, 11h30, 11h45
#ECH = PSCOPF.generate_ech(network, TS, mode) #ech: -4h, -1h, -30mins, -15mins, 0h

#sequence = PSCOPF.generate_sequence(network, TS, ECH, mode)

# Personalised sequence

# tout est ferme
sequence_ferme = PSCOPF.Sequence(Dict([
    ts1 - Dates.Minute(45)  => [PSCOPF.BalanceMarket(), PSCOPF.TSOBilevel()],
]))

# Décisions d'impositions fermes
sequence_impositions_ferme = PSCOPF.Sequence(Dict([
    ts1 - Dates.Minute(45)  => [PSCOPF.BalanceMarket(), PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL=true))],
    ts1 - Dates.Minute(15)  => [PSCOPF.Assessment()],
]))

# tout est par scénario
sequence_free = PSCOPF.Sequence(Dict([
    ts1 - Dates.Minute(45)  => [PSCOPF.BalanceMarket(), PSCOPF.TSOBilevel()],
    ts1 - Dates.Minute(15)  => [PSCOPF.Assessment()],
]))


PSCOPF.rm_non_prefixed(output_path, "pscopf_")
exec_context = PSCOPF.PSCOPFContext(network, TS, mode,
                                    generators_init_state,
                                    uncertainties, nothing,
                                    output_path)

PSCOPF.run!(exec_context, sequence_impositions_ferme)
