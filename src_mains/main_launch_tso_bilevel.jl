"""
    main_launch_tso_bilevel

Parameters:
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

matpower_case = length(ARGS) > 0 ? ARGS[1] : "case5"
output_folder = joinpath(@__DIR__, "..", "data", matpower_case)

###############################################
# CONFIGS :
###############################################
PSCOPF.set_config!("CONSIDER_N_1",false)
# Activate dynamic solving
PSCOPF.set_config!("ADD_RSO_CSTR_DYNAMICALLY", true)
# LOL sensibility problems
PSCOPF.set_config!("tso_loss_of_load_penalty_value", 1e4)
PSCOPF.set_config!("market_loss_of_load_penalty_value", 1e4)
PSCOPF.set_config!("big_m_value", 2.e4)
PSCOPF.set_config!("lol_eps", 0.01)
PSCOPF.set_config!("FIX_NULL_LOL", true)
# LOG options
PSCOPF.set_config!("LOG_COMBINATIONS", true)
PSCOPF.set_config!("LOG_NB_CSTRS_FILENAME", joinpath(@__DIR__, "..", "main_nbCSTRS.log"))
PSCOPF.set_config!("TEMP_GLOBAL_LOGFILE", joinpath(@__DIR__, "..", "main_timingsRSO.log"))
PSCOPF.init_logging(joinpath(output_folder, "main_launch_market_tsobilevel.log"))


###############################################
# INPUT & PARAMS
###############################################

ts1 = DateTime("2015-01-01T11:00:00")
ECH = [ts1-Hour(4), ts1-Hour(2), ts1-Hour(1), ts1-Minute(30), ts1-Minute(15), ts1]
TS = PSCOPF.create_target_timepoints(ts1)

#################################################################################################################
# Launch
#################################################################################################################

###############################################
# Read instance
###############################################
instance_path = joinpath(output_folder, "instance")
generated_network = PSCOPF.Data.pscopfdata2network(instance_path)
uncertainties = PSCOPF.PSCOPFio.read_uncertainties(instance_path)
gen_init = PSCOPF.PSCOPFio.read_initial_state(instance_path)

logfile = PSCOPF.get_config("TEMP_GLOBAL_LOGFILE")
open(logfile, "a") do file_l
    write(file_l, "-"^120 * "\n")
    write(file_l, @sprintf("usecase : %s\n", output_folder))
    write(file_l, @sprintf("dynamic? : %s\n", PSCOPF.get_config("ADD_RSO_CSTR_DYNAMICALLY")))
    write(file_l, @sprintf("n-1? : %s\n", PSCOPF.get_config("CONSIDER_N_1")))
    write(file_l, @sprintf("nb rso constraints : %d\n", PSCOPF.nb_rso_constraint(generated_network, length(PSCOPF.get_scenarios(uncertainties)), length(TS))))
end


###############################################
# Solve usecase : mode 1
###############################################
output_path = joinpath(output_folder, "tso_bilevel_30Min")
mode = PSCOPF.ManagementMode("tsobilevel", Minute(60))
sequence = PSCOPF.Sequence(Dict([
        ts1 - Dates.Minute(4*60)  => [PSCOPF.EnergyMarket()],
        ts1 - Dates.Minute(30)  => [PSCOPF.TSOBilevel()],
        ts1 - Dates.Minute(15)  => [PSCOPF.Assessment()],
    ]))

PSCOPF.rm_non_prefixed(output_path, "pscopf_")
exec_context = PSCOPF.PSCOPFContext(generated_network, TS, mode,
                                    gen_init,
                                    uncertainties, nothing,
                                    output_path)
PSCOPF.run!(exec_context, sequence)

println(PSCOPF.TIMER_TRACKS)
