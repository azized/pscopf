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

matpower_case = length(ARGS) > 0 ? ARGS[1] : "case4gs"

function main_launch_tso_bilevel(case_p::String;
            dynamic_p::Bool = true, nb_cstr_add_per_iter_p::Int = 10, n_1_p::Bool = true, null_lol_p::Bool = true)
    TimerOutputs.reset_timer!(PSCOPF.TIMER_TRACKS)

    output_folder = joinpath(@__DIR__, "..", "data", case_p)


    ###############################################
    # CONFIGS :
    ###############################################
    PSCOPF.set_config!("CONSIDER_N_1",n_1_p)
    # Activate dynamic solving
    PSCOPF.set_config!("ADD_RSO_CSTR_DYNAMICALLY", dynamic_p)
    PSCOPF.set_config!("MAX_ADD_RSO_CSTR_PER_ITER", nb_cstr_add_per_iter_p)
    # LOL sensibility problems
    PSCOPF.set_config!("FIX_NULL_LOL", null_lol_p)
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

end

main_launch_tso_bilevel(matpower_case)
