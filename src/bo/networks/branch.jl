using DataStructures

Limits = SortedDict{String,Float64} # network_case => limit

struct Branch
    id::String
    #src::Int # bus id from which emits the branch
    #dst::Int # bus it to which goes the branch

    # Metier
    limit::Limits
end

function Branch(id::String, limit::Float64)
    limits = Limits(Networks.BASECASE=>limit)
    return Branch(id, limits)
end

function get_limit(branch::Branch)
    return branch.limit
end

function safeget_limit(branch::Branch, network_case::String)::Float64
    limits = get_limit(branch)
    if !haskey(limits, network_case)
        msg_l = @sprintf("missing limit value for branch %s in the network case %s", get_id(branch), network_case)
        error(msg_l)
    end
    return limits[network_case]
end

function add_limit!(branch::Branch, network_case::String, limit::Float64)
    branch.limit[network_case] = limit
end

################
## INFO / LOG ##
################

function get_info(branch::Branch)::String
    info::String =
        string(branch.id)
        # * ":" *
        # string(branch.src) *
        # "->" *
        # string(branch.dst)
    return info
end
