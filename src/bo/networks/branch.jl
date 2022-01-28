
struct Branch
    src::Int # bus id from which emits the branch
    dst::Int # bus it to which goes the branch
    # Metier
end


################
## INFO / LOG ##
################

function get_info(branch::Branch)::String
    info::String = 
        string(branch.src) *
        "->" *
        string(branch.dst)
    return info
end
