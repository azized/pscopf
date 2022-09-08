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
# include(joinpath(root_path, "src", "PTDF.jl"));

function main_launch_tso_bilevel(case_p::String;
            dynamic_p::Bool = true,
            nb_cstr_add_per_iter_p::Int = 10, viols_add_fct::String = "TS",
            n_1_p::Bool = true, null_lol_p::Bool = true
            )
    TimerOutputs.reset_timer!(PSCOPF.TIMER_TRACKS)

    input_folder = joinpath(@__DIR__, "..", "data", case_p)
    name_l = "tso_bilevel_"*string(dynamic_p)*"_"*viols_add_fct
    output_path = joinpath(input_folder, name_l)


    ###############################################
    # CONFIGS :
    ###############################################
    PSCOPF.set_config!("CONSIDER_N_1",n_1_p)
    # Activate dynamic solving
    PSCOPF.set_config!("ADD_RSO_CSTR_DYNAMICALLY", dynamic_p)
    PSCOPF.set_config!("MAX_ADD_RSO_CSTR_PER_ITER", nb_cstr_add_per_iter_p)
    PSCOPF.set_config!("DYNAMIC_ONLY_STEP1", true)
    PSCOPF.set_config!("VIOLATIONS_ADD_FUNCTION", viols_add_fct)
    # LOL sensibility problems
    PSCOPF.set_config!("FIX_NULL_LOL", null_lol_p)
    # LOG options
    PSCOPF.set_config!("EXTRA_LOG", true)
    PSCOPF.set_config!("LOG_COMBINATIONS", true)
    PSCOPF.set_config!("LOG_NB_CSTRS_FILENAME", joinpath(@__DIR__, "..", "main_nbCSTRS.log"))
    PSCOPF.set_config!("TEMP_GLOBAL_LOGFILE", joinpath(@__DIR__, "..", "main_timingsRSO.log"))
    PSCOPF.init_logging(joinpath(output_path, "main_launch_market_tsobilevel.log"))


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
    instance_path = joinpath(input_folder, "instance")
    generated_network = PSCOPF.Data.pscopfdata2network(instance_path)
    uncertainties = PSCOPF.PSCOPFio.read_uncertainties(instance_path)
    gen_init = PSCOPF.PSCOPFio.read_initial_state(instance_path)

    logfile = PSCOPF.get_config("TEMP_GLOBAL_LOGFILE")
    open(logfile, "a") do file_l
        write(file_l, "-"^120 * "\n")
        write(file_l, @sprintf("usecase : %s\n", output_path))
        write(file_l, @sprintf("dynamic? : %s\n", PSCOPF.get_config("ADD_RSO_CSTR_DYNAMICALLY")))
        write(file_l, @sprintf("n-1? : %s\n", PSCOPF.get_config("CONSIDER_N_1")))
        write(file_l, @sprintf("nb rso constraints : %d\n", PSCOPF.nb_rso_constraint(generated_network, length(PSCOPF.get_scenarios(uncertainties)), length(TS))))
    end


    ###############################################
    # Solve usecase : mode 1
    ###############################################
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

function upon_termination()
    @info PSCOPF.TIMER_TRACKS
end

function str_to_bool(str)
    if lowercase(str) == "true"
        return true
    elseif lowercase(str) == "false"
        return false
    else
        error("cannot convert string to bool : "*str)
    end
end

atexit(upon_termination)

matpower_case = length(ARGS) > 0 ? ARGS[1] : "case89pegase"
dynamic = length(ARGS) > 1 ? str_to_bool(ARGS[2]) : true
add_fct = length(ARGS) > 2 ? ARGS[3] : "TS"
max_add_per_iter = (add_fct in ["TS", "TS_GROUP"]) ? 2 : 8

main_launch_tso_bilevel(matpower_case,
                        dynamic_p=dynamic, nb_cstr_add_per_iter_p=max_add_per_iter, viols_add_fct=add_fct)
