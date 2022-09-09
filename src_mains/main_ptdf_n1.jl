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


MATPOWER_NETWORKS = [
    "case4gs",
    "case5",
    "case6ww",
    "case9",
    "case9Q",
    "case9target",
    "case14",
    "case24_ieee_rts",
    "case30",
    "case30pwl",
    "case30Q",
    "case39",
    "case57",
    "case89pegase",
    "case118",
    "case145",
    "case300",
    # "case1354pegase",
    # "case13659pegase",
    # "case1888rte",
    # "case1951rte",
    # "case2383wp",
    # "case2736sp",
    # "case2737sop",
    # "case2746wop",
    # "case2746wp",
    # "case2848rte",
    # "case2868rte",
    # "case2869pegase",
    # "case3012wp",
    # "case3120sp",
    # "case3375wp",
    # "case6468rte",
    # "case6470rte",
    # "case6495rte",
    # "case6515rte",
    # "case9241pegase",
];

#########################
# EXECUTION
#########################
function main(input_path, ref_bus_num, distributed, eps_diag)
    network = PTDF.read_network(input_path)
    out_path = input_path
    PTDF.compute_and_write_n_non_bridges(network, ref_bus_num, distributed, eps_diag, input_path, out_path)
end

input_path = ( length(ARGS) > 0 ? ARGS[1] :
                    joinpath(@__DIR__, "..", "data", "ptdf", "3buses_3branches") )
ref_bus_num = 1
distributed = true
eps_diag = 1e-3

IN = length(ARGS) > 0 ? [ARGS[1]] : MATPOWER_NETWORKS
for matpower_case in IN
    @info matpower_case
    matpower_path = joinpath(@__DIR__, "..", "data_matpower", matpower_case)
    main(matpower_path, ref_bus_num, distributed, eps_diag)
end
