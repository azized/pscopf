module Networks

# using ..AmplTxt

include("./generator.jl")
include("./bus.jl")
include("./branch.jl")
include("./network.jl")

function get_id(obj)
    return obj.id
end

export
    # struct
    Generator, Bus, Branch, Network,
    PTDFDict, PTDFValues, Limits,
    # enum
    GeneratorType,
    # functions
    ## bus
    add_new_bus!, add_new_buses!, add_bus!, add_buses!,
    get_bus, safeget_bus, get_buses, get_nb_buses,
    get_generators,
    ## branch
    add_new_branch!, add_new_branches!, add_branch!, add_branches!,
    add_new_limit!,
    get_branch, safeget_branch, get_branches,
    get_limit, safeget_limit,
    ## load and generator
    add_new_generator_to_bus!,
    get_generator, safeget_generator, get_generators, get_generators_of_type,
    get_bus_id, get_type, get_p_min, get_p_max, get_start_cost, get_prop_cost, get_dmo, get_dp,
    needs_commitment, PILOTABLE, LIMITABLE,
    ## ptdf
    add_ptdf_elt!,
    get_ptdf_component, safeget_ptdf_component,
    get_ptdf, safeget_ptdf,
    ## utils
    get_id,
    safeget_generator_or_bus, get_generator_or_bus,
    get_info

end #module Networks
