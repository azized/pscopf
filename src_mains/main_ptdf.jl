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

MATPOWER_NETWORKS = [
    "case118"
    ,"case1354pegase"
    ,"case13659pegase"
    ,"case14"
    ,"case145"
    ,"case1888rte"
    ,"case1951rte"
    ,"case2383wp"
    ,"case24_ieee_rts"
    ,"case2736sp"
    ,"case2737sop"
    ,"case2746wop"
    ,"case2746wp"
    ,"case2848rte"
    ,"case2868rte"
    ,"case2869pegase"
    ,"case30"
    ,"case300"
    ,"case3012wp"
    ,"case30pwl"
    ,"case30Q"
    ,"case3120sp"
    ,"case3375wp"
    ,"case39"
    ,"case4gs"
    ,"case5"
    ,"case57"
    ,"case6468rte"
    ,"case6470rte"
    ,"case6495rte"
    ,"case6515rte"
    ,"case6ww"
    ,"case89pegase"
    ,"case9"
    ,"case9241pegase"
    ,"case9Q"
    ,"case9target"
];

function compute_ptdf(input_path::String, eps_diag_p=1e-6)
    output_path = joinpath(input_path, "pscopf_ptdf.txt")
    ref_bus_num = 1
    distributed = true

    #########################
    # EXECUTION
    #########################
    network = PTDF.read_network(input_path)
    ptdf = PTDF.compute_ptdf(network, ref_bus_num, eps_diag_p)
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

# input_path = ( length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "..", "data_matpower", "case1354pegase") )
for input_path in MATPOWER_NETWORKS
    compute_ptdf(joinpath(@__DIR__, "..", "data_matpower", input_path));
end
