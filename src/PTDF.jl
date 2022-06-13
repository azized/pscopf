#TODO tests

module PTDF

using DataStructures
using Parameters
using Printf
using SparseArrays, LinearAlgebra

include("utils.jl")

mutable struct Bus
    id::Int
    name::String
end

mutable struct Branch
    from::Int
    to::Int
    r::Float64
    x::Float64
    name::String
end
function get_b(branch::Branch)
    # return branch.x / sqrt(branch.r * branch.r + branch.x * branch.x)
    return 1/branch.x;
end

@with_kw mutable struct Network
    bus_to_i::Dict{Int,Int} = Dict{Int,Int}()
    branch_to_i::Dict{Int,Int} = Dict{Int,Int}()

    branches::Dict{Int,Branch} = Dict{Int,Branch}()
    buses::Dict{Int,Bus} = Dict{Int,Bus}()
    CC0::Dict{Int,Int} = Dict{Int,Int}()
end

function read_network(dir_path)
    network = Network();
    read_buses!(network, dir_path);
    read_branches!(network, dir_path);

    return network
end

function get_or_add!(d::Dict{Int,Int}, i::Int)
    if ! haskey(d, i)
        next_id = length(d) + 1;
        d[i] = next_id;
    end
    return d[i];
end

function read_buses!(network::Network, dir_path)
    open(joinpath(dir_path, "buses.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = split_with_space(ln);

                name = buffer[1];
                num = parse(Int, buffer[2]);
                numCC = parse(Int, buffer[3]);
                if numCC == 0
                    id = get_or_add!(network.bus_to_i, num);
                    push!(network.buses, id => Bus(num, name));
                    network.CC0[num] = 1;
                end
            end
        end
    end
end


function cut_branch!(network::Network, i_branch_to_cut::Int)
    cut_branch = pop!(network.branches, i_branch_to_cut)
    @info("removed branch ", cut_branch.name, " from network!")
    #no deletion in branch_to_i to keep the idx and to show the branch in the output but with zeros
    #change branch_to_i ?
    return cut_branch, network
end

function reduced_network(network::Network, i_branch_to_cut::Int)
    # if a new CC is generated even the EOD constraint is no longer global,
    #the two induced networks function independantly
    @warn("TODO : verify that the network is still connected and no new CC was generated!")
    network_l = deepcopy(network)

    return cut_branch!(network_l, i_branch_to_cut)
end



function read_branches!(network::Network, dir_path)
    open(joinpath(dir_path, "branches.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = split_with_space(ln);

                name = buffer[1]
                num  = parse(Int, buffer[2]);
                num_or = parse(Int, buffer[3]);
                num_ex = parse(Int, buffer[4]);

                on_CC0 = haskey(network.CC0, num_or) && haskey(network.CC0, num_ex);

                if on_CC0 && num_or != -1  && num_ex != -1
                    id = get_or_add!(network.branch_to_i, num);
                    ior = get_or_add!(network.bus_to_i, num_or);
                    iex = get_or_add!(network.bus_to_i, num_ex);
                    r = parse(Float64, buffer[5]);
                    x = parse(Float64, buffer[6]);
                    push!(network.branches, id => Branch(ior, iex, r, x, name));
                end
            end
        end
    end
end


function get_i_for_bus(network_p::Network, bus_name_p::String)
    bus_id_l = first(filter( x -> x[2].name == bus_name_p , network_p.buses))[2].id
    println(bus_name_p, " ; id ", bus_id_l, " ; i ", network_p.bus_to_i[bus_id_l])
    println("BUSES:\n", network_p.buses)
    return network_p.bus_to_i[bus_id_l]
    #or for now, return first(filter( x -> x[2].name == bus_name_l , network.buses))[1]
end

function get_B(network::Network, DIAG_EPS::Float64)
    B_dict = Dict{Tuple{Int,Int},Float64}();
    for kvp in collect(network.branches)
        b = get_b(kvp[2]);

        ior = kvp[2].from;
        iex = kvp[2].to;


        # ior->iex
        # + dans ior
        # - dans iex

        # diagonal part
        # sur ior
        # +theta_ior
        # -theta_iex
        B_dict[ior, ior] = get(B_dict, (ior, ior), DIAG_EPS) + b;
        B_dict[ior, iex] = get(B_dict, (ior, iex), 0) - b;

        # sur iex
        # -theta_ior
        # +theta_iex
        B_dict[iex, ior] = get(B_dict, (iex, ior), 0) - b;
        B_dict[iex, iex] = get(B_dict, (iex, iex), DIAG_EPS) + b;


        # B_dict[ior, ior] = get(B_dict, (ior, ior), DIAG_EPS) + b;
        # B_dict[iex, iex] = get(B_dict, (iex, iex), DIAG_EPS) - b;

        # B_dict[ior, iex] = get(B_dict, (ior, iex), 0) + b;
        # B_dict[iex, ior] = get(B_dict, (iex, ior), 0) - b;
    end
    Brow = Int[];
    Bcol = Int[];
    Bval = Float64[];
    for kvp in collect(B_dict)
        push!(Brow, kvp[1][1]);
        push!(Bcol, kvp[1][2]);
        push!(Bval, kvp[2]);
    end
    n = length(network.bus_to_i);
    println("n = ", n)
    B = sparse(Brow, Bcol, Bval, n, n);
    Bdense = Matrix(B);
    return Bdense
end

function write_B(file_path_p, B::Matrix, network_p::Network)
    nb_buses_l = size(B)[1]
    @assert( size(B)[1] == size(B)[2] == length(network_p.bus_to_i) )

    open(file_path_p, "w") do file_l
        write(file_l, @sprintf("#%20s %20s %20s\n", "BUS_1", "BUS_2", "b"))
        for i_bus1_l in 1:nb_buses_l, i_bus2_l in 1:nb_buses_l
            bus1_name_l =  @sprintf("\"%s\"", network.buses[i_bus1_l].name)
            bus2_name_l =  @sprintf("\"%s\"", network.buses[i_bus2_l].name)
            write(file_l, @sprintf("%20s %20s %20.6E\n", bus1_name_l, bus2_name_l, B[i_bus1_l, i_bus2_l]))
        end
    end
end

function get_B_inv(Bdense::Matrix, ref_bus::Int64)
    println("ref_bus is ", ref_bus)
    n = size(Bdense)[1]
    Bdense2 = Bdense
    Bdense2[ref_bus, 1:n] .= 0
    Bdense2[1:n, ref_bus] .= 0
    Bdense2[ref_bus, ref_bus] = 1
    binv = inv(Bdense2);
    binv[ref_bus,ref_bus] = 0;
    return binv
end

function get_PTDF(network::Network, binv::Matrix, ref_bus::Int)
    n = length(network.bus_to_i);
    m = length(network.branch_to_i);
    println("n ", n)
    println("m ", m)
    PTDF = zeros(Float64, m, n);
    for kvp in collect(network.branches)
        b = get_b(kvp[2]);

        branchid = kvp[1];
        ior = kvp[2].from;
        iex = kvp[2].to;

        for i in 1:n
            PTDF[branchid, i] += b * binv[ior, i];
            PTDF[branchid, i] -= b * binv[iex, i];
            # if ior != ref_bus
            #     PTDF[branchid, i] += b * binv[ior, i];
            # end
            # if iex != ref_bus
            #     PTDF[branchid, i] -= b * binv[iex, i];
            # end
        end
    end
    println("PTDF=\n",PTDF)
    return PTDF
end

function distribute_slack(PTDF::Matrix, coeffs::Vector{Float64})
    m, n = size(PTDF);
    v = PTDF * coeffs;
    result = -v*(ones(n)');
    result +=PTDF ;
    return result;
end

function distribute_slack(PTDF::Matrix, coeffs_p::Dict{String,Float64}, network_p::Network)
    m, n = size(PTDF);
    @assert( n == length(coeffs_p) )

    vect_coeffs_l = zeros(n)
    for (bus_name_l, coeff_l) in coeffs_p
        i_l = get_i_for_bus(network_p, bus_name_l)
        vect_coeffs_l[i_l] = coeff_l
    end

    #normalize coeffs
    vect_coeffs_l = vect_coeffs_l / sum(vect_coeffs_l)

    return distribute_slack(PTDF, vect_coeffs_l)
end

function distribute_slack(PTDF::Matrix)
    m, n = size(PTDF);
    μ = 1/n;
    coeffs_l = ones(n)*μ
    return distribute_slack(PTDF, coeffs_l)
end

function write_slack_distribution(file_path::String, network::Network, coeffs_p::Vector{Float64})
    n = length(network.bus_to_i);
    m = length(network.branches);
    open(file_path, "w") do file
        write(file, @sprintf("#%20s %20s\n", "BUS", "COEFF"))
        for bus_i_l in 1:n
            bus_name =  @sprintf("\"%s\"", network.buses[bus_i_l].name)
            write(file, @sprintf("%20s %20.6E\n", bus_name, coeffs_p[bus_i_l]))
        end
    end
end

function write_PTDF(file_path::String,
                    network::Network, PTDF::Matrix,
                    distributed=false, ref_bus::Int=-1,
                    i_cut_branch::Union{Nothing,Int}=nothing;
                    PTDF_TRIMMER::Float64=1e-06,)
    n = length(network.bus_to_i);
    m = length(network.branch_to_i);
    open(file_path, "w") do file
        ref_name =  distributed ? "\"distributed\"" : @sprintf("\"%s\"", network.buses[ref_bus].name)
        cut_branch =  isnothing(i_cut_branch) ? "\"NONE\"" : @sprintf("\"%s\"", network.branches[i_cut_branch].name)
        write(file, @sprintf("#%20s %20s %20s %20s\n", "REF_BUS", ref_name, "CUT_BRANCH", cut_branch))
        for branch_id in 1:m
            for bus_id in 1:n
                branch_name =  @sprintf("\"%s\"", network.branches[branch_id].name)
                bus_name =  @sprintf("\"%s\"", network.buses[bus_id].name)
                ptdf = PTDF[branch_id,bus_id]
                if abs(ptdf)<PTDF_TRIMMER
                    ptdf = 0.
                end
                write(file, @sprintf("%20s %20s %20.6E\n", branch_name, bus_name,ptdf))
            end
        end
    end
end

function compute_ptdf(network, ref_bus::Int=1)
    B = get_B(network, 1e-6);
    Binv = get_B_inv(B, ref_bus);
    PTDF = get_PTDF(network, Binv, ref_bus);

    return PTDF
end

export compute_ptdf

export get_b;
export get_B;
export get_B_inv;
export get_PTDF;
export write_PTDF;
export distribute_slack;


end #module PTDF
