using PSCOPF

using Test
using Dates
using DataStructures

@testset verbose=true "example_usecase" begin

    #=
    S: [S1,S2]
    TS: [11h, 11h30]
                      FO
          EM   EM   EM |
    ECH: [7h, 7h30, 10h, ...]
                            bus 1                   bus 2
                               |                      |
                               |                      |
           TS1         TS2     |                      |       TS1         TS2
                               |                      |
        load-------------------|                      |----load
    S1:   ,   ,   |   ,   ,    |                      |S1:   ,   ,   |   ,   ,
    S2:   ,   ,   |   ,   ,    |                      |S2:   ,   ,   |   ,   ,
                               |                      |
                               |                      |
        wind_1-----------------|                      |----wind_2
             0, 200            |                      |         0, 200
             0, 10             |                      |         0, 11
        15mins, 15mins         |                      |    15mins, 15mins
    S1:   ,   ,   |   ,   ,    |                      |S1:   ,   ,   |   ,   ,
    S2:   ,   ,   |   ,   ,    |                      |S2:   ,   ,   |   ,   ,
                               |                      |
                               |                      |
        ccg_1------------------|                      |----ccg_2
           150, 600            |                      |       100, 600
           450, 10             |                      |       50k, 20
            4h, 15mins         |                      |        4h, 15mins
                               |                      |
                               |                      |
        tac_1------------------|                      |
            10, 300            |                      |
           120, 20             |                      |
        30mins, 15mins         |                      |
                               |                      |

    =#
    @testset "example_small_usecase_with_reserves_units" begin
        out_path = joinpath(@__DIR__, "..", "..", "default_out", "example_small_usecase")
        rm(out_path, recursive=true, force=true)

        ECH = [DateTime("2015-01-01T07:00:00"), DateTime("2015-01-01T07:30:00"),DateTime("2015-01-01T10:00:00"),DateTime("2015-01-01T11:00:00")]
        TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:30:00")]

        network = PSCOPF.Networks.Network()
        # Buses
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        PSCOPF.Networks.add_new_bus!(network, "bus_2")
        # Branches
        PSCOPF.Networks.add_new_branch!(network, "branch_1_2", 500.);
        # PTDF
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_1", 0.5)
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_2", -0.5)
        #Generators - Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1", PSCOPF.Networks.LIMITABLE,
                                                0., 200.,
                                                0., 1.,
                                                Dates.Second(15*60), Dates.Second(15*60)) #dmo, dp : always 0. ?
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "wind_2", PSCOPF.Networks.LIMITABLE,
                                                0., 200.,
                                                0., 2.,
                                                Dates.Second(15*60), Dates.Second(15*60))
        #Generators - Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "ccg_1", PSCOPF.Networks.PILOTABLE,
                                                150., 600., #pmin, pmax
                                                450., 10., #start_cost, prop_cost
                                                Dates.Second(4*3600), Dates.Second(15*60)) #dmo, dp
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "tac_1", PSCOPF.Networks.PILOTABLE,
                                                10., 300.,
                                                120., 20.,
                                                Dates.Second(30*60), Dates.Second(15*60))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "ccg_2", PSCOPF.Networks.PILOTABLE,
                                                100., 600., #pmin, pmax
                                                500., 20., #start_cost, prop_cost
                                                Dates.Second(4*3600), Dates.Second(15*60)) #dmo, dp
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "reserve_1", PSCOPF.Networks.PILOTABLE,
                                                0., 600., #pmin, pmax
                                                0., 500., #start_cost, prop_cost
                                                Dates.Second(0*3600), Dates.Second(0*60)) #dmo, dp
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "reserve_2", PSCOPF.Networks.PILOTABLE,
                                                0., 600., #pmin, pmax
                                                0., 500., #start_cost, prop_cost
                                                Dates.Second(0*3600), Dates.Second(0*60)) #dmo, dp

        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1", TS[1], "S1", 610.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1", TS[2], "S1", 590.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_2", TS[1], "S1", 110.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_2", TS[2], "S1", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1", TS[1], "S1", 49.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1", TS[2], "S1", 50.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_2", TS[1], "S1", 78.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_2", TS[2], "S1", 80.)
        #S2
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1", TS[1], "S2", 640.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1", TS[2], "S2", 615.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_2", TS[1], "S2", 110.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_2", TS[2], "S2", 105.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1", TS[1], "S2", 18.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1", TS[2], "S2", 23.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_2", TS[1], "S2", 62.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_2", TS[2], "S2", 58.)
        # Load increased, wind prod decreased => start tac_1 (not firm)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_1", TS[1], "S1", 630.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_1", TS[2], "S1", 600.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_2", TS[1], "S1", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_2", TS[2], "S1", 110.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_1", TS[1], "S1", 40.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_1", TS[2], "S1", 41.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_2", TS[1], "S1", 81.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_2", TS[2], "S1", 78.)
        #S2
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_1", TS[1], "S2", 620.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_1", TS[2], "S2", 610.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_2", TS[1], "S2", 120.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_2", TS[2], "S2", 110.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_1", TS[1], "S2", 35.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_1", TS[2], "S2", 31.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_2", TS[1], "S2", 60.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_2", TS[2], "S2", 62.)
        # Finally the increase is not that high => don't start tac_1
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_1", TS[1], "S1", 615.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_1", TS[2], "S1", 613.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_2", TS[1], "S1", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_2", TS[2], "S1", 97.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_1", TS[1], "S1", 42.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_1", TS[2], "S1", 40.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_2", TS[1], "S1", 78.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_2", TS[2], "S1", 75.)
        #S2
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_1", TS[1], "S2", 620.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_1", TS[2], "S2", 610.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_2", TS[1], "S2", 95.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "bus_2", TS[2], "S2", 90.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_1", TS[1], "S2", 53.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_1", TS[2], "S2", 50.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_2", TS[1], "S2", 76.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[3], "wind_2", TS[2], "S2", 75.)


        # initial generators state
        generators_init_state = SortedDict(
                        "ccg_1" => PSCOPF.OFF,
                        "ccg_2" => PSCOPF.OFF,
                        "tac_1" => PSCOPF.OFF,
                    )
        mode = PSCOPF.PSCOPF_MODE_1
        sequence = PSCOPF.generate_sequence(network, TS, ECH, mode)
        exec_context = PSCOPF.PSCOPFContext(network, TS, mode, generators_init_state,
                                            uncertainties, nothing, out_path)
        PSCOPF.run!(exec_context, sequence)
    end

    println("\n\n\n\n\n")


    @testset "example_2bus" begin
        network = PSCOPF.Networks.Network()

        # Buses
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        PSCOPF.Networks.add_new_bus!(network, "bus_2")

        # Branches
        PSCOPF.Networks.add_new_branch!(network, "branch_1_2", 500.);

        # PTDF
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_1", 0.5)
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_2", -0.5)
        #Alternatively,
        # network.ptdf = SortedDict{String,SortedDict{String, Float64}}(
        #                         "branch_1_2" => SortedDict{String, Float64}(
        #                             "bus_1" => 0.5,
        #                             "bus_2" => -0.5,
        #                             ),
        #                         )

        # Generators
        #Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_0", PSCOPF.Networks.LIMITABLE,
                                                0., 0., #pmin, pmax : Not concerned ? min is always 0, max is the limitation
                                                0., 1., #start_cost, prop_cost : start cost is always 0 ?
                                                Dates.Second(0), Dates.Second(0)) #dmo, dp : always 0. ?
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "wind_2_0", PSCOPF.Networks.LIMITABLE,
                                                0., 0.,
                                                0., 2.,
                                                Dates.Second(0), Dates.Second(0))
        #Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "ccg_1_0", PSCOPF.Networks.PILOTABLE,
                                                10., 200., #pmin, pmax
                                                450., 10., #start_cost, prop_cost
                                                Dates.Second(4*3600), Dates.Second(15*60)) #dmo, dp
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "tac_2_0", PSCOPF.Networks.PILOTABLE,
                                                10., 200.,
                                                120., 20.,
                                                Dates.Second(30*60), Dates.Second(15*60))
        # initial generators state
        generators_init_state = SortedDict(
                        "ccg_1_0" => PSCOPF.ON,
                        "tac_2_0" => PSCOPF.OFF,
                    )

        # Uncertainties
        # uncertainties = PSCOPF.PSCOPFio.read_uncertainties(data_folder) #filename must be : pscopf_uncertainties.txt
        uncertainties = PSCOPF.Uncertainties(
                                DateTime("2015-01-01T07:00:00") => SortedDict(
                                        "wind_1_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "wind_2_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_1" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_2" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                ),
                                DateTime("2015-01-01T10:00:00") => SortedDict(
                                        "wind_1_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "wind_2_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_1" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_2" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                ),
                                DateTime("2015-01-01T10:30:00") => SortedDict(
                                        "wind_1_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "wind_2_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_1" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_2" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                ),
                                DateTime("2015-01-01T10:45:00") => SortedDict(
                                        "wind_1_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "wind_2_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_1" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_2" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                ),
                                DateTime("2015-01-01T11:00:00") => SortedDict(
                                        "wind_1_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "wind_2_0" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_1" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                        "bus_2" => SortedDict(
                                            DateTime("2015-01-01T11:00:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:15:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:30:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                            DateTime("2015-01-01T11:45:00") => SortedDict(
                                                                                        "S1" => 0.,
                                                                                        "S2" => 0.,
                                                                                        ),
                                        ),
                                ),
        )
        #Alternatively, PSCOPF.add_uncertainty!(uncertainties, ech, nodal_injection_name, ts, scenario_name, value)

        # Timesteps
        TS = [DateTime("2015-01-01T11:00:00"),
                DateTime("2015-01-01T11:15:00"),
                DateTime("2015-01-01T11:30:00"),
                DateTime("2015-01-01T11:45:00")]

        mode = PSCOPF.PSCOPF_MODE_1
        ECH = PSCOPF.generate_ech(network, TS, mode)

        sequence = PSCOPF.generate_sequence(network, TS, ECH, mode)

        exec_context = PSCOPF.PSCOPFContext(network, TS, mode, generators_init_state, uncertainties, nothing)
        PSCOPF.run!(exec_context, sequence)

        market_schedule = PSCOPF.get_market_schedule(exec_context)
        tso_schedule = PSCOPF.get_tso_schedule(exec_context)
    end

end
