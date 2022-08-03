using PSCOPF

using Test
using Dates
using DataStructures

@testset verbose=true "test_flows" begin

    #=
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
    Csta=0, Cprop=2     |                      | Csta=0, Cprop=3
    ON                  |                      | ON
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10     S1: 15 |                      | S1: 40  S1: 50
      S2: 10     S2: 15 |                      | S2: 45  S2: 50

    =#

    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T07:00:00")
    network = PSCOPF.Networks.Network()
    # Buses
    PSCOPF.Networks.add_new_bus!(network, "bus_1")
    PSCOPF.Networks.add_new_bus!(network, "bus_2")
    # Branches
    PSCOPF.Networks.add_new_branch!(network, "branch_1_2", 35.);
    # PTDF
    PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_1", 0.5)
    PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_2", -0.5)
    # Limitables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                            0., 100.,
                                            0., 1.,
                                            Dates.Second(0), Dates.Second(0))
    # Pilotables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                            0., 100.,
                                            0., 2.,
                                            Dates.Second(0), Dates.Second(0))
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                            0., 100.,
                                            0., 3.,
                                            Dates.Second(0), Dates.Second(0))
    # Uncertainties
    uncertainties = PSCOPF.Uncertainties()
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 20.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S2", 30.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:15:00"), "S1", 15.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:15:00"), "S2", 30.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 10.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 10.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S1", 15.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 15.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", DateTime("2015-01-01T11:00:00"), "S1", 40.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", DateTime("2015-01-01T11:00:00"), "S2", 45.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", DateTime("2015-01-01T11:15:00"), "S1", 50.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", DateTime("2015-01-01T11:15:00"), "S2", 50.)

    function create_schedule()
        schedule = PSCOPF.Schedule(PSCOPF.Utilitary(), DateTime("2015-01-01T06:00:00"), SortedDict(
            "wind_1_1" => PSCOPF.GeneratorSchedule("wind_1_1",
                SortedDict(),
                SortedDict(TS[1]=> PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>20.,"S2"=>30.)),
                            TS[2] => PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>15.,"S2"=>30.)),
                            ),
                ),
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(),
                SortedDict(TS[1]=> PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>30.,"S2"=>25.)),
                            TS[2] => PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>50.,"S2"=>35.)),
                            ),
                ),
            "prod_2_1" => PSCOPF.GeneratorSchedule("prod_2_1",
                SortedDict(),
                SortedDict(TS[1]=> PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>0.,"S2"=>0.)),
                            TS[2] => PSCOPF.UncertainValue{Float64}(missing, SortedDict("S1"=>0.,"S2"=>0.)),
                            ),
                ),
            )
        )

        return  schedule
    end

    @testset "test_flows_simple" begin
        schedule = create_schedule()
        @test 40. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S1")
        @test 45. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S2")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S1")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S2")
    end

    @testset "test_flows_simple_change" begin
        schedule = create_schedule()
        PSCOPF.set_prod_value!(schedule, "prod_1_1", TS[1], "S1", 10.) # 30 -> 10 : -20
        PSCOPF.set_prod_value!(schedule, "prod_2_1", TS[1], "S1", 20.) # 0 -> 20 : +20

        # changed
        @test 20. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S1")

        # These are unchanged
        @test 45. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S2")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S1")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S2")
    end

    @testset "test_flows_slack" begin
        schedule = create_schedule()
        PSCOPF.set_prod_value!(schedule, "prod_1_1", TS[1], "S1", 10.) # 30 -> 10 : -20
        PSCOPF.set_prod_value!(schedule, "prod_2_1", TS[1], "S1", 10.) # 0 -> 10 : +10

        # slack effect
        @test ( PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S1")
                ≈ ( 0.5 * ( (20. + 10.) - (10.) ) #bus1
                   -0.5 * ( (10.) - (40.) ) #bus2
                )#25
            )

        # These are unchanged
        @test 45. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S2")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S1")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S2")
    end

    @testset "test_flows_capping" begin
        schedule = create_schedule()
        PSCOPF.set_capping_value!(schedule, "wind_1_1", TS[1], "S1", 5.)

        # capping 5MW from wind_1_1 => 20-5
        @test ( PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S1")
                ≈ ( 0.5 * ( ((20. - 5.) + 30.) - (10.) ) #bus1
                   -0.5 * ( (0.) - (40.) ) #bus2
                )#37.5
            )

        # These are unchanged
        @test 45. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S2")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S1")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S2")
    end

    @testset "test_flows_loss_of_load" begin
        schedule = create_schedule()
        PSCOPF.set_loss_of_load_value!(schedule, "bus_2", TS[1], "S1", 15.)

        # cut a 15MW conso on bus2 => 40-15
        @test ( PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S1")
                ≈ ( 0.5 * ( (20. + 30.) - (10.) ) #bus1
                   -0.5 * ( (0.) - (40. - 15.) ) #bus2
                )#32.5
            )

        # These are unchanged
        @test 45. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[1], "S2")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S1")
        @test 50. ≈ PSCOPF.compute_flow("branch_1_2", uncertainties, schedule, network, ech, TS[2], "S2")
    end

end
