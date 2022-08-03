using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_tso_constraints" begin

    #=
        ECH = 10h
        TS: [11h, 11h15]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 20    S1: 15    |----------------------|
        S2: 30    S2: 30    |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
    INIT: ON                |                      |INIT: ON
                            |                      | DP=DMO=2h
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 10     S1: 17   |                      | S1: 40  S1: 48
        S2: 10     S2: 13   |                      | S2: 45  S2: 52
    =#

    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T07:00:00")

    tso = PSCOPF.TSOOutFO()
    context = PSCOPFFixtures.context_2buses_2TS_2S(TS, ech)

    next_ech = DateTime("2015-01-01T07:30:00") # all decisions will be free
    result, firmness = PSCOPF.run_step!(context, tso, ech, next_ech)

    @testset "tso_successful_launch" begin
        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
    end

    @testset "tso_updated_schedule_respects_the_input_firmness" begin
        @test context.tso_schedule.decision_time == ech
        @test PSCOPF.verify_firmness(firmness, context.tso_schedule,
                        excluded_ids=PSCOPF.get_limitables_ids(context))
    end

    @testset "tso_respects_EOD" begin
        @test context.tso_schedule.decision_time == ech
        @test PSCOPF.verify_firmness(firmness, context.tso_schedule,
                        excluded_ids=PSCOPF.get_limitables_ids(context))
        for ts in TS
            for s in ["S1", "S2"]
                @test ( abs( PSCOPF.compute_eod(PSCOPF.get_uncertainties(context),
                                        PSCOPF.get_tso_schedule(context),
                                        PSCOPF.get_network(context),
                                        ech, ts, s) )
                        < 1e-09 )
            end
        end
    end

    @testset "tso_respects_RSO_constraints" begin
        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        #prod_2_1 is used for RSO constraints
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_2_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_2_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_2_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_2_1", TS[2], "S2")
        #prod_1_1 might be used to deviate the least from market schedule

        # Note : compute_flow does not consider loss_of_load (valid here cause loss_of_load=0, proven by pscopf_OPTIMAL)
        for ts in TS
            for s in ["S1", "S2"]
                flow = PSCOPF.compute_flow("branch_1_2",
                                        PSCOPF.get_uncertainties(context),
                                        PSCOPF.get_tso_schedule(context),
                                        PSCOPF.get_network(context),
                                        ech, ts, s)
                @test PSCOPF.in_bounds(flow, -35., 35.)
            end
        end
    end

end
