"""
    main_ptdf
A main file to generate a PTDF file

Parameters:
    input_path : path to input data directory describing a grid
                (not the pscopf_ files but branches.txt and buses.txt)
"""

using Dates
using Printf

root_path = dirname(@__DIR__)
push!(LOAD_PATH, root_path);
cd(root_path)
include(joinpath(root_path, "src", "PTDF.jl"));




#########################
# INPUT & PARAMS
#########################

input_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "data", "ptdf", "3buses_3branches") )
ref_bus_num = 1
distributed = true


function compute(network_p, ref_bus_num_p, distributed_p, i_cut_branch_p=nothing)
    if isnothing(i_cut_branch_p)
        network_l = network_p
        name = ""
    else
        cut_branch_l, network_l = PTDF.reduced_network(network_p, i_cut_branch_p)
        name = "_n1_"*cut_branch_l.name
    end

    ptdf_l = PTDF.compute_ptdf(network_l, ref_bus_num_p)
    if distributed_p
        ptdf_l = PTDF.distribute_slack(ptdf_l);
        # coeffs = Dict([ "poste_1_0" => .2,
        #                 "poste_2_0" => .8])
        # ptdf_l = PTDF.distribute_slack(ptdf_l, coeffs, network_l);
    end
    filename = @sprintf("pscopf_ptdf%s.txt", name)
    output_path = joinpath(input_path, filename)
    #use original network to print full ptdf containing the cut branch (with 0 coeffs)
    PTDF.write_PTDF(output_path, network_p, ptdf_l, distributed_p, ref_bus_num_p, i_cut_branch_p)
end

#########################
# EXECUTION
#########################
network = PTDF.read_network(input_path)

compute(network, ref_bus_num, distributed)
for (i_cut_branch,_) in network.branch_to_i
    compute(network, ref_bus_num, distributed, i_cut_branch)
end
