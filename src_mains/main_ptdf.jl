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




#########################
# INPUT & PARAMS
#########################

input_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "data", "ptdf", "2buses") )
output_path = joinpath(input_path, "pscopf_ptdf.txt")
ref_bus_num = 1
distributed = true

#########################
# EXECUTION
#########################
network = PTDF.read_network(input_path)
ptdf = PTDF.compute_ptdf(network, ref_bus_num)
if distributed
    ptdf = PTDF.distribute_slack(ptdf);
    # coeffs = Dict([ "poste_1_0" => .2,
    #                 "poste_2_0" => .8])
    # ptdf = PTDF.distribute_slack(ptdf, coeffs, network);
end
PTDF.write_PTDF(output_path, network, ptdf, distributed, ref_bus_num)
