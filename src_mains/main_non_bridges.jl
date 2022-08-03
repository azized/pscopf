"""
    main_non_bridges
A main file to check which branches are non-bridges
i.e. branches whom removal does not result in two non-connected grids

Parameters:
    input_path : path to input data directory describing a grid
                (not the pscopf_ files but branches.txt and buses.txt)
"""

using Graphs
using MetaGraphs
using Printf

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

function network2graph(network::PTDF.Network)::Tuple{MetaGraph, Set{PTDF.Branch}}
    graph::MetaGraph = MetaGraph(SimpleGraph())
    parallel_lines::Set{PTDF.Branch} = Set()

    for i in range(1, length(network.buses))
        bus = network.buses[i]
        # vertex i is bus i
        if !add_vertex!(graph, :bus_ref, bus)
            msg_l = @sprintf("error adding vertex of bus %s to the graph", bus.name)
            error(msg_l)
        end
    end

    for (_, branch) in network.branches
        if has_edge(graph, branch.from, branch.to) # Undirected graph
            push!(parallel_lines, branch)

            existing_edge = first(filter_edges(graph,
                                    (g,e) -> ((src(e)==branch.from && dst(e)==branch.to) || (src(e)==branch.to && dst(e)==branch.from)) ))
            existing_branch = get_prop(graph, existing_edge, :branch_ref)
            push!(parallel_lines, existing_branch)
        else
            add_edge!(graph, branch.from, branch.to, :branch_ref, branch)
        end
    end

    return graph, parallel_lines
end

function is_multiline(edge, graph, parallel_lines)::Bool
    return (get_prop(graph, edge, :branch_ref) in parallel_lines)
end

function write_non_bridges(graph::MetaGraph, parallel_lines::Set{PTDF.Branch}, file_path::String)
    bridges_l = bridges(graph)
    # only keep non parallel lines (cause parallel lines cannot be bridges)
    filter!(e -> !is_multiline(e, graph, parallel_lines), bridges_l)
    @info @sprintf("graph contains %d bridges!", length(bridges_l))


    non_bridges::Set{PTDF.Branch} = Set(parallel_lines)
    for e in edges(graph)
        if !(e in bridges_l)
            push!(non_bridges, get_prop(graph, e, :branch_ref))
        end
    end

    mkpath(dirname(file_path))
    open(file_path, "w") do file
        write(file, @sprintf("#NON_BRIDGES  #nb_bridges=%d\n", length(bridges_l)) )
        for branch in non_bridges
            write(file, branch.name*"\n")
        end
    end
end


function compute_non_bridges(input_path::String)
    network = PTDF.read_network(input_path)

    graph, parallel_lines = network2graph(network)

    output_path = joinpath(input_path, "non_bridges.txt")
    write_non_bridges(graph, parallel_lines, output_path)
end


#########################
# MAIN
#########################
function main(dir_names)
    for input_dirname in dir_names
        @info input_dirname
        input_path = joinpath(@__DIR__, "..", "data_matpower", input_dirname)
        compute_non_bridges(input_path);
    end
end


# MATPOWER_NETWORKS = ["case14"]
main(MATPOWER_NETWORKS)

# input_path = joinpath(@__DIR__, "..", "data", "ptdf", "3buses_3branches")
# compute_non_bridges(input_path);



