module FlowCalculator

PTDF_FILE = "ptdf.txt"
INJECTIONS_FILE = "injections.txt"
CONSUMPTIONS_FILE = "consumptions.txt"

using Parameters
using Printf
using DataStructures

struct FC_Injection
    bus::String
    value::Float64
end

struct FC_Consumption
    bus::String
    value::Float64
end

@with_kw struct FC_Bus
    name::String=""
    injections::SortedDict{String,FC_Injection}=SortedDict{String,FC_Injection}()
    consumptions::SortedDict{String,FC_Consumption}=SortedDict{String,FC_Consumption}()
end
function FC_Bus(name::String)
    return FC_Bus(name, SortedDict{String,FC_Injection}(), SortedDict{String,FC_Consumption}())
end

struct FC_Branch
    name::String
end

#usecase,Branch,Bus_id => value)
FC_PTDF = Dict{Tuple{String,String,String},Float64}

@with_kw struct FC_Network
    buses::SortedDict{String,FC_Bus}=SortedDict{String,FC_Bus}()
    branches::SortedDict{String,FC_Branch}=SortedDict{String,FC_Branch}()
    ptdf::FC_PTDF=FC_PTDF()
    network_cases::Set{String}=Set{String}()
end

function split_with_space(str::String)
    result = String[];
    if length(str) > 0
        start_with_quote = startswith(str, "\"");
        buffer_quote = split(str, keepempty=false, "\"");
        i = 1;
        while i <= length(buffer_quote)
            if i > 1 || !start_with_quote
                str2 = buffer_quote[i];
                buffer_space = split(str2, keepempty=false);
                for str3 in buffer_space
                    push!(result, str3);
                end
                i += 1;
            end
            if i <= length(buffer_quote)
                push!(result, buffer_quote[i]);
                i += 1;
            end
        end
    end
    return result;
end

function read_network(data)
    network = FC_Network()

    #read ptdf
    open(joinpath(data, PTDF_FILE), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = split_with_space(ln);
                ptdf_case = buffer[1]
                branch_id = buffer[2]
                bus_id = buffer[3]
                ptdf_value = parse(Float64, buffer[4])
                network.ptdf[ptdf_case, branch_id, bus_id] = ptdf_value
                push!(network.network_cases, ptdf_case)
                get!(network.branches, branch_id, FC_Branch(branch_id))
            end
        end
    end

    #read consumptions
    open(joinpath(data, CONSUMPTIONS_FILE), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = split_with_space(ln);
                conso_name = buffer[1]
                bus_name = buffer[2]
                conso_value = parse(Float64, buffer[3])

                bus_l = get!(network.buses, bus_name, FC_Bus(bus_name))
                if haskey(bus_l.consumptions, conso_name)
                    error("duplicate consumption name")
                end
                bus_l.consumptions[conso_name] = FC_Consumption(bus_name, conso_value)
            end
        end
    end

    #read injections
    open(joinpath(data, INJECTIONS_FILE), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = split_with_space(ln);
                inj_name = buffer[1]
                bus_name = buffer[2]
                inj_value = parse(Float64, buffer[3])

                bus_l = get!(network.buses, bus_name, FC_Bus(bus_name))
                if haskey(bus_l.injections, inj_name)
                    error("duplicate injection name")
                end
                bus_l.injections[inj_name] = FC_Injection(bus_name, inj_value)
            end
        end
    end

    return network
end

function compute_flow(network::FC_Network,
                    branch_name::String, ptdf_case::String)
    flow_l = 0.
    for (bus_name,bus) in network.buses
        ptdf_val = network.ptdf[ptdf_case, branch_name, bus_name]

        # + injections
        for (_, inj) in bus.injections
            flow_l += ptdf_val * inj.value
        end

        # - loads
        for (_, conso) in bus.consumptions
            flow_l -= ptdf_val * conso.value
        end
    end

    return flow_l
end

function compute_flows(network::FC_Network)
    @printf("%20s %20s %s\n", "case", "branch", "flow")
    for case in network.network_cases
        for (_,branch) in network.branches
            flow_l = compute_flow(network, branch.name, case)
            @printf("%20s %20s %f\n", case, branch.name, flow_l)
        end
    end
end

function compute_flows(dir::String)
    network = read_network(dir)
    compute_flows(network)
end

export compute_flows

end # module FlowCalculator


in_data = joinpath(@__DIR__, "..", "data", "flowcalculator")
FlowCalculator.compute_flows(in_data)
