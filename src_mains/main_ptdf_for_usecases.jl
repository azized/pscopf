"""
    main_ptdf
A main file to generate a PTDF file

Parameters:
    input_path : path to input data directory describing a grid
                (not the pscopf_ files but branches.txt and buses.txt)
"""

using Dates

root_path = dirname(@__DIR__)
push!(LOAD_PATH, root_path);
cd(root_path)
include(joinpath(root_path, "src", "PTDF.jl"));


function compute_ptdf(input_path::String, ref_bus_num::Int=1, EPS_DIAG=1e-6, distributed=false)    
    output_path = joinpath(input_path, "pscopf_ptdf.txt")

    #########################
    # EXECUTION
    #########################
    network = PTDF.read_network(input_path)
    ptdf = PTDF.compute_ptdf(network, ref_bus_num, EPS_DIAG)
    if distributed
        ptdf = PTDF.distribute_slack(ptdf);
        # coeffs = Dict([ "poste_1_0" => .2,
        #                 "poste_2_0" => .8])
        # ptdf = PTDF.distribute_slack(ptdf, coeffs, network);
    end
    PTDF.write_PTDF(output_path, network, ptdf, distributed, ref_bus_num)

end


#########################
# INPUT & PARAMS
#########################

compute_ptdf(joinpath(@__DIR__, "..", "usecases-euro-simple", "usecase2-tunnel-puissance", "data"), 3)