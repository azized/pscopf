using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_energy_market_at_fo" begin

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

    prod_2_1's level was already decided for TS1 and TS2 : 10 and 15MW

    =#


    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T10:00:00")
    context = PSCOPFFixtures.context_2buses_2TS_2S(TS, ech)
    context.tso_schedule = PSCOPF.Schedule(PSCOPF.Market(), Dates.DateTime("2015-01-01T09:00:00"), SortedDict(
                                            "wind_1_1" => PSCOPF.GeneratorSchedule("wind_1_1",
                                                SortedDict{Dates.DateTime, PSCOPF.UncertainValue{PSCOPF.GeneratorState}}(),
                                                SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                                                        SortedDict("S1"=>20., "S2"=>30.)),
                                                            TS[2] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                                                        SortedDict("S1"=>15., "S2"=>30.))),
                                                ),
                                            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                                                SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                        SortedDict("S1"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                                                            TS[2] => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                        SortedDict("S1"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                                                            ),
                                                SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                                                        SortedDict("S1"=>20., "S2"=>15.)),
                                                            TS[2] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                                                        SortedDict("S1"=>35., "S2"=>20.))
                                                            )
                                                ),
                                            "prod_2_1" => PSCOPF.GeneratorSchedule("prod_2_1",
                                                SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                        SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                            TS[2] => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                        SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                            ),
                                                SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(10.,
                                                                                                                        SortedDict("S1"=>10., "S2"=>10.)),
                                                            TS[2] => PSCOPF.UncertainValue{PSCOPF.Float64}(15.,
                                                                                                                        SortedDict("S1"=>15., "S2"=>15.))
                                                            )
                                                )
                                            )
                                    )

    market = PSCOPF.EnergyMarketAtFO()
    next_ech = DateTime("2015-01-01T10:01:00") # does not matter : all decisions will be DECIDED or TODECIDE
    result, firmness = PSCOPF.run_step!(context, market, ech, next_ech)

    @testset "energy_market_at_fo_is_launched_at_fo" begin
        @test ech == (TS[1] - PSCOPF.get_fo_length(PSCOPF.get_management_mode(context)))
    end

    @testset "energy_market_required_decisions_are_all_firm" begin
        println(firmness)
        expected_firmness = PSCOPF.Firmness(
                        SortedDict("prod_1_1" => SortedDict(TS[1] => PSCOPF.TO_DECIDE,
                                                            TS[2] => PSCOPF.TO_DECIDE,
                                                            ),
                                    "prod_2_1" => SortedDict(TS[1] => PSCOPF.DECIDED,
                                                            TS[2] => PSCOPF.DECIDED,
                                                            ),
                                    ),
                        SortedDict("wind_1_1" => SortedDict(TS[1] => PSCOPF.TO_DECIDE,
                                                            TS[2] => PSCOPF.TO_DECIDE,
                                                            ),
                                    "prod_1_1" => SortedDict(TS[1] => PSCOPF.TO_DECIDE,
                                                            TS[2] => PSCOPF.TO_DECIDE,
                                                            ),
                                    "prod_2_1" => SortedDict(TS[1] => PSCOPF.DECIDED,
                                                            TS[2] => PSCOPF.DECIDED,
                                                            ),
                                    )
                        )
        @test expected_firmness == firmness
    end

    @testset "energy_market_at_fo_successful_launch" begin
        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
    end

    #=
    EnergyMarketAtFO returns firm decisions
    =#
    @testset "energy_market_at_fo_all_decisions_are_firm" begin
        @test context.market_schedule.decision_time == ech
        for (gen_id,gen_schedule) in context.market_schedule.generator_schedules
            if PSCOPF.Networks.get_type(PSCOPF.Networks.get_generator(PSCOPF.get_network(context), gen_id)) == PSCOPF.Networks.PILOTABLE
                for (ts, uncertain) in gen_schedule.production
                    @test PSCOPF.is_definitive(uncertain)
                end
            end
        end
        for (gen_id,gen_schedule) in context.market_schedule.generator_schedules
            for (ts, uncertain) in gen_schedule.commitment
                @test PSCOPF.is_definitive(uncertain)
            end
        end
    end

    #=
    EnergyMarketAtFO operates on an aggregated scenario
    =#
    @testset "energy_market_at_fo_aggregated_scenario" begin
        @test PSCOPF.SCENARIOS_DELIMITER == "_+_"

        expected_s_agg = "S1_+_S2"
        s_agg = PSCOPF.aggregate_scenario_name(context, ech)
        @test s_agg == expected_s_agg
    end

    #=
    EnergyMarketAtFO operates on an aggregated scenario and respects preceding decisions
    preceding decision : prod_2_1 started at 10MW for TS[1]

    The aggregate scenario :
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
        25      22.5    |----------------------|
                        |         35           |
                        |                      |
                        |                      |
    (pilotable) prod_1_1|                      |(pilotable) prod_2_1
    Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
    Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
    ON    ON       ON   |                      | ON  ON      ON
                        |                      |
           load(bus_1)  |                      |load(bus_2)
        10      15      |                      |    42.5  50

    =#
    @testset "energy_market_at_fo_production_levels" begin
        s_agg = PSCOPF.aggregate_scenario_name(context, ech)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        @test 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S2")
        @test 17.5 ≈ value(result.pilotable_model.p_injected["prod_1_1", TS[1], s_agg])
        @test 10. ≈ value(result.pilotable_model.p_injected["prod_2_1", TS[1], s_agg])
        @test value(result.limitable_model.p_global_capping[TS[1], s_agg]) < 1e-09
        @test value(result.lol_model.p_global_loss_of_load[TS[1], s_agg]) < 1e-09

        @test 22.5 ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[2], "S1")
        @test 22.5 ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[2], "S2")
        @test 27.5 ≈ value(result.pilotable_model.p_injected["prod_1_1", TS[2], s_agg])
        @test 15. ≈ value(result.pilotable_model.p_injected["prod_2_1", TS[2], s_agg])
        @test value(result.limitable_model.p_global_capping[TS[2], s_agg]) < 1e-09
        @test value(result.lol_model.p_global_loss_of_load[TS[2], s_agg]) < 1e-09

        @testset "energy_market_at_fo_respect_preceding_decisions" begin
            #prod_2_1 was already set to 10 for TS1
            @test PSCOPF.ON == PSCOPF.safeget_commitment_value(context.market_schedule, "prod_2_1", TS[1])
            @test 10. ≈ PSCOPF.safeget_prod_value(context.market_schedule, "prod_2_1", TS[1])
            #prod_2_1 was already set to 15 for TS2
            @test PSCOPF.ON == PSCOPF.safeget_commitment_value(context.market_schedule, "prod_2_1", TS[2])
            @test 15. ≈ PSCOPF.safeget_prod_value(context.market_schedule, "prod_2_1", TS[2])
        end
    end

    #=
    EnergyMarketAtFO considers aggregated scenario uncertainties,
    => It does not necessarily respect the EOD of the initial individual scenarios
    If we look at the initial scenarios we can have : overloads, overproductions or balanced scenarios
    =#
    @testset "energy_market_at_fo_does_not_necessarily_respect_EOD_on_initial_scenarios" begin
        s_agg = PSCOPF.aggregate_scenario_name(context, ech)

        #TS[1], "S1" : aggregated scenarios induces an over prodction
        @test PSCOPF.is_different(0. ,PSCOPF.compute_eod(PSCOPF.get_uncertainties(context),
                                                        PSCOPF.get_market_schedule(context),
                                                        PSCOPF.get_network(context),
                                                        ech, TS[1], "S1"))
        #TS[1], "S2" : aggregated scenarios induces an over load
        @test PSCOPF.is_different(0., PSCOPF.compute_eod(PSCOPF.get_uncertainties(context),
                                                        PSCOPF.get_market_schedule(context),
                                                        PSCOPF.get_network(context),
                                                        ech, TS[1], "S2") )
        #TS[2] : initial scenarios are balanced
        @test PSCOPF.compute_eod(PSCOPF.get_uncertainties(context),
                                PSCOPF.get_market_schedule(context),
                                PSCOPF.get_network(context),
                                ech, TS[2], "S1") < 1e-09
        @test PSCOPF.compute_eod(PSCOPF.get_uncertainties(context),
                                PSCOPF.get_market_schedule(context),
                                PSCOPF.get_network(context),
                                ech, TS[2], "S2") < 1e-09

        # capping and LoL are defined for the aggregated scenario
        @test value( result.lol_model.p_global_loss_of_load[TS[1], s_agg] ) < 1e-09
        @test value( result.limitable_model.p_global_capping[TS[1], s_agg] ) < 1e-09
        @test value( result.lol_model.p_global_loss_of_load[TS[2], s_agg] ) < 1e-09
        @test value( result.limitable_model.p_global_capping[TS[2], s_agg] ) < 1e-09
    end

    #=
    EnergyMarketAtFO satisfies the EOD of the aggregated scenario
    NOTE: If the EOD could not be satisfied, the market might use capping and LoL for the aggregated scenario
    =#
    @testset "energy_market_at_fo_respects_EOD_on_aggregated_uncertainties" begin
        s_agg, agg_uncertainties = PSCOPF.aggregate_scenarios(context, ech)
        agg_uncertainties = PSCOPF.Uncertainties(ech => agg_uncertainties)
        for ts in TS
            @test_skip 1e-09 > abs( PSCOPF.compute_eod(agg_uncertainties,
                                                PSCOPF.get_market_schedule(context),
                                                PSCOPF.get_network(context),
                                                ech, ts, s_agg) )
            @test value( result.lol_model.p_global_loss_of_load[ts, s_agg] ) < 1e-09
            @test value( result.limitable_model.p_global_capping[ts, s_agg] ) < 1e-09
        end
    end

end
