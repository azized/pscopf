module PSCOPFFixtures

using PSCOPF

using Dates
using DataStructures
using Printf

"""
creates a network with 2 buses
"""
function network_2buses(;limit::Float64=35., ptdf1=0.5, ptdf2=-0.5)
    network = PSCOPF.Networks.Network()
    # Buses
    PSCOPF.Networks.add_new_bus!(network, "bus_1")
    PSCOPF.Networks.add_new_bus!(network, "bus_2")
    # Branches
    PSCOPF.Networks.add_new_branch!(network, "branch_1_2", limit);
    # PTDF
    PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_1", ptdf1)
    PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_2", ptdf2)

    return network
end

#=
    ECH = 10h
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 20    S1: 15  |----------------------|
      S2: 30    S2: 30  |         35           |
                        |                      |
                        |                      |
    (pilotable) prod_1_1|                      |(pilotable) prod_2_1
    Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
    Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
INIT: ON                |                      |INIT: ON
                        |                      | DP=DMO=2h
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10     S1: 17 |                      | S1: 40  S1: 48
      S2: 10     S2: 13 |                      | S2: 45  S2: 52
=#
"""
Creates a context with preset uncertainties and for a preset 2 buses network :
    - on bus_1 : a limitable and an pilotable
    - on bus_2 : an pilotable
    pilotable unit of bus_2 has a DP=DMO=2h
Arguments :
    TS  : a 2 elements DateTime vector for the target timesteps
    ech : DateTime giving the current horizon timepoint
"""
function context_2buses_2TS_2S(TS, ech;
                                logs=nothing)
    @assert(length(TS) == 2)
    network = PSCOPFFixtures.network_2buses()

    # Limitables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                            0., 100.,
                                            0., 1.,
                                            Dates.Second(0), Dates.Second(0))
    # Pilotables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            450., 10.,
                                            Dates.Second(0), Dates.Second(0))
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            800., 15.,
                                            Dates.Second(2*60*60), Dates.Second(2*60*60))
    # initial generators state
    generators_init_state = SortedDict(
                    "prod_1_1" => PSCOPF.ON,
                    "prod_2_1" => PSCOPF.ON
                )
    # Uncertainties
    uncertainties = PSCOPF.Uncertainties()
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], "S1", 20.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], "S2", 30.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[2], "S1", 15.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[2], "S2", 30.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", 10.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S2", 10.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[2], "S1", 17.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[2], "S2", 13.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", 40.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S2", 45.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[2], "S1", 48.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[2], "S2", 52.)

    context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                generators_init_state,
                                uncertainties, nothing,
                                logs)

    return context
end



end #module PSCOPFFixtures

