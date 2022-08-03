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

function press_to_continue()
    disable_timer!(PSCOPF.TIMER_TRACKS)
    println("\n"^3)
    println("press enter to continue")
    readline()
    println("\n"^3)
    enable_timer!(PSCOPF.TIMER_TRACKS)
end


###############################################
# INPUT & PARAMS
###############################################

matpower_case = length(ARGS) > 0 ? ARGS[1] : "case5"
input_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "data_matpower", matpower_case) )
output_folder = joinpath(@__DIR__, "..", "data", matpower_case)

ts1 = DateTime("2015-01-01T11:00:00")
ECH = [ts1-Hour(4), ts1-Hour(2), ts1-Hour(1), ts1-Minute(30), ts1-Minute(15), ts1]
TS = PSCOPF.create_target_timepoints(ts1)

#################################################################################################################
# Launch
#################################################################################################################

###############################################
# CONFIGS :
###############################################
# Activate dynamic solving
PSCOPF.set_config!("ADD_RSO_CSTR_DYNAMICALLY", true)
# LOL sensibility problems
PSCOPF.set_config!("tso_loss_of_load_penalty_value", 1e6)
PSCOPF.set_config!("market_loss_of_load_penalty_value", 1e6)
PSCOPF.set_config!("big_m_value", 2.e6)
PSCOPF.set_config!("lol_eps", 0.01)
# LOG options
PSCOPF.set_config!("LOG_COMBINATIONS", true)
PSCOPF.set_config!("LOG_NB_CSTRS_FILENAME", joinpath(@__DIR__, "..", "main_nbCSTRS.log"))
PSCOPF.set_config!("TEMP_GLOBAL_LOGFILE", joinpath(@__DIR__, "..", "main_timingsRSO.log"))
PSCOPF.init_logging(joinpath(output_folder, "main_launch_1_2.log"))


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
output_mode_1 = joinpath(output_folder, "mode1")
mode_1 = PSCOPF.PSCOPF_MODE_1
sequence = PSCOPF.generate_sequence(generated_network, TS, ECH, mode_1)

PSCOPF.rm_non_prefixed(output_mode_1, "pscopf_")
exec_context = PSCOPF.PSCOPFContext(generated_network, TS, mode_1,
                                    gen_init,
                                    uncertainties, nothing,
                                    output_mode_1)
time_mode_1 = @elapsed begin
    try
        PSCOPF.run!(exec_context, sequence)
    catch e
        showerror(stdout, e)
    end
end

# press_to_continue()
# open(logfile, "a") do file_l
#     write(file_l, "-"^60 * "\n")
# end


###############################################
# Solve usecase : mode 2
###############################################
output_mode_2 = joinpath(output_folder, "mode2")
mode_2 = PSCOPF.PSCOPF_MODE_2
sequence = PSCOPF.generate_sequence(generated_network, TS, ECH, mode_2)

PSCOPF.rm_non_prefixed(output_mode_2, "pscopf_")
exec_context = PSCOPF.PSCOPFContext(generated_network, TS, mode_2,
                                    gen_init,
                                    uncertainties, nothing,
                                    output_mode_2)
time_mode_2 = @elapsed begin
    try
        PSCOPF.run!(exec_context, sequence)
    catch e
        showerror(stdout, e)
    end
end



println("Mode 1 took :", time_mode_1)
println("Mode 2 took :", time_mode_2)

println(PSCOPF.TIMER_TRACKS)


open(logfile, "a") do file_l
    write(file_l, @sprintf("Mode 1 took : %s\n", time_mode_1))
    write(file_l, @sprintf("Mode 2 took : %s\n", time_mode_2))
    write(file_l, @sprintf("TIMES:\n%s\n", PSCOPF.TIMER_TRACKS))
end
