using PSCOPF

using Test
using Dates
using DataStructures

@testset verbose=true "test_schedule" begin

    @testset "test_market_empty_schedule" begin
        println("\n\n\n")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus");
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus", "gen_1", PSCOPF.Networks.IMPOSABLE, 0., 0., 0., 0., Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus", "gen_2", PSCOPF.Networks.LIMITABLE, 0., 0., 0., 0., Dates.Second(0), Dates.Second(0))

        TS = [DateTime("2015-01-01T11:00:00"),
                DateTime("2015-01-01T11:15:00"),
                DateTime("2015-01-01T11:30:00"),
                DateTime("2015-01-01T11:45:00")]
        scenarios = ["S1","S2"]
        ech = DateTime("2015-01-01T07:00:00")
        schedule = PSCOPF.Schedule(PSCOPF.Market(), ech)
        PSCOPF.init!(schedule, network, TS, scenarios)

        @test PSCOPF.is_market(schedule.decider)
        @test !PSCOPF.is_tso(schedule.decider)
        @test schedule.decision_time == ech
        @test length(schedule.values) == 4
        for ts in TS
            @test length(PSCOPF.get_values(schedule, ts)) == 2
            for generator in PSCOPF.Networks.get_generators(network)
                @test ismissing(PSCOPF.get_value(schedule, ts, PSCOPF.Networks.get_id(generator)))
            end
        end
    end

end
