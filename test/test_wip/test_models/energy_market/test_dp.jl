using PSCOPF

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_energy_market_dp" begin

    #=
    TS: [11h]
    S: [S1]
                        bus 1
                        |
    (limitable) wind_1_1|load
    Pmin=0, Pmax=100    | S1: 55
    Csta=0, Cprop=1     |
    DP => 10h40         |
      S1: 20            |
                        |
    (pilotable) prod_1_1|
    Pmin=10, Pmax=100   |
    Csta=0, Cprop=10    |
    DP => 10h30         |
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=0, Cprop=15   |
     DP => 10h45        |
                        |

    =#
    function create_instance(ech, next_ech, ts,
                            wind_1_level,
                            prod_1_level,
                            prod_2_level;
                            fo_length=Minute(5))
        TS = [ts]
        ECH = [DateTime("2015-01-01T10:00:00"), DateTime("2015-01-01T10:40:00")]
        network = PSCOPF.Networks.Network()
        # Buses
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(20*60), Dates.Second(20*60))
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                10., 100.,
                                                0., 10.,
                                                Dates.Second(3*60*60), Dates.Second(30*60))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_2", PSCOPF.Networks.PILOTABLE,
                                                10., 100.,
                                                0., 15.,
                                                Dates.Second(3*60*60), Dates.Second(15*60))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1_1", TS[1], "S1", 20.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1", TS[1], "S1", 55.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "wind_1_1", TS[1], "S1", 20.)
        PSCOPF.add_uncertainty!(uncertainties, ECH[2], "bus_1", TS[1], "S1", 55.)
        # initial generators state : need to pay starting cost at TS[1]
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        #ManagementMode
        mode = PSCOPF.ManagementMode("test_mode", fo_length)

        context = PSCOPF.PSCOPFContext(network, TS, mode,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        firmness = PSCOPF.compute_firmness(market,
                                            ech, next_ech,
                                            TS, context)

        prod_1_definitive_state = (prod_1_level < 1e-09) ? PSCOPF.OFF : PSCOPF.ON
        prod_1_definitive_level = (PSCOPF.get_power_level_firmness(firmness, "prod_1_1", ts) == PSCOPF.DECIDED) ? prod_1_level : missing
        prod_2_definitive_state = (prod_2_level < 1e-09) ? PSCOPF.OFF : PSCOPF.ON
        prod_2_definitive_level = (PSCOPF.get_power_level_firmness(firmness, "prod_1_2", ts) == PSCOPF.DECIDED) ? prod_2_level : missing
        context.tso_schedule = PSCOPF.Schedule(PSCOPF.Market(), Dates.DateTime("2015-01-01T06:00:00"), SortedDict(
                "wind_1_1" => PSCOPF.GeneratorSchedule("wind_1_1",
                    SortedDict{Dates.DateTime, PSCOPF.UncertainValue{PSCOPF.GeneratorState}}(),
                    # SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                    #                                                                         SortedDict("S1"=>PSCOPF.ON))),
                    SortedDict(ts => PSCOPF.UncertainValue{Float64}(missing,
                                                                    SortedDict("S1"=>wind_1_level)))
                    ),
                "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                    SortedDict(ts => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(prod_1_definitive_state,
                                                                    SortedDict("S1"=>prod_1_definitive_state))),
                    SortedDict(ts => PSCOPF.UncertainValue{Float64}(prod_1_definitive_level,
                                                                    SortedDict("S1"=>prod_1_level)))
                    ),
                "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                    SortedDict(ts => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(prod_2_definitive_state,
                                                                    SortedDict("S1"=>prod_2_definitive_state))),
                    SortedDict(ts => PSCOPF.UncertainValue{Float64}(prod_2_definitive_level,
                                                                    SortedDict("S1"=>prod_2_level)))
                    )
                )
            )

        return context, market, firmness
    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1
                        |
    (limitable) wind_1_1|load
    Pmin=0, Pmax=100    | S1: 55
    Csta=0, Cprop=1     |
    DP => 10h40         |
      S1: 20            |
                        |
    (pilotable) prod_1_1|
    Pmin=10, Pmax=100   |
    Csta=0, Cprop=10    |
    DP => 10h30         |
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=0, Cprop=15   |
     DP => 10h45        |
                        |

    ECH1*         ECH2            ECH3                                 <---FO--->TS
    |             |               |
    10h           10h30           10h40          10h45                          11h
                  <------------------------------------------DP(prod1)----------->
                                  <--------------------------DP(wind1)----------->
                                                   <---------DP(prod2)----------->

    We can change the production level of a unit before reaching the DP
    ech = 10h
    TS = 11h
    DP(prod_1_1) = 30 mins => 10h30
    DP(prod_1_2) = 15 mins => 10h45

    We suppose that in the preceding step we decided :
    prod_1_1 : 0. (OFF)
    prod_1_2 : 37. (ON)

    In terms of cost, it is cheaper to use prod_1_1. But, since we are past DMO, we can no longer start prod_1_1.
    wind_1_1 : 20 MW
    remaining Demand : 35MW
    Since we haven't reached the DP yet, we can change the production level of prod_1_2 from 37 MW to 35 MW
    => the new decision :
    wind_1_1 : 20 MW
    prod_1_1 : OFF (unchanged)
    prod_1_2 : 35. (decreased from 37MW)
    =#
    @testset "energy_market_can_decrease_the_production_level_before_dp_if_ON" begin
        ech = DateTime("2015-01-01T10:00:00")
        next_ech = DateTime("2015-01-01T10:30:00")
        TS = [DateTime("2015-01-01T11:00:00")]
        context, market, firmness = create_instance(ech, next_ech, TS[1],
                                                    18., 0., 37.,
                                                    )
        expected_firmness = PSCOPF.Firmness(
                SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.DECIDED),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.DECIDED), ),
                SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                            "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE), )
                )
        @test firmness == expected_firmness

        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        # prod_1_1 is OFF cause it already was OFF:
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        # prod_1_2 was ON, it can change it's production : 37 => 35.
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1
                        |
    (limitable) wind_1_1|load
    Pmin=0, Pmax=100    | S1: 55
    Csta=0, Cprop=1     |
    DP => 10h40         |
      S1: 20            |
                        |
    (pilotable) prod_1_1|
    Pmin=10, Pmax=100   |
    Csta=0, Cprop=10    |
    DP => 10h30         |
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=0, Cprop=15   |
     DP => 10h45        |
                        |

    ECH1*         ECH2            ECH3                                 <---FO--->TS
    |             |               |
    10h           10h30           10h40          10h45                          11h
                  <------------------------------------------DP(prod1)----------->
                                  <--------------------------DP(wind1)----------->
                                                   <---------DP(prod2)----------->

    We can change the production level of a unit before reaching the DP
    ech = 10h
    TS = 11h
    DP(prod_1_1) = 30 mins => 10h30
    DP(prod_1_2) = 15 mins => 10h45

    We suppose that in the preceding step we decided :
    prod_1_1 : 0. (OFF)
    prod_1_2 : 33. (ON)

    In terms of cost, it is cheaper to use prod_1_1. But, since we are past DMO, we can no longer start prod_1_1.
    wind_1_1 : 20 MW
    remaining Demand : 35MW
    Since we haven't reached the DP yet, we can change the production level of prod_1_2 from 33 MW to 35 MW
    => the new decision :
    wind_1_1 : 20 MW
    prod_1_1 : OFF (unchanged)
    prod_1_2 : 35. (increased from 33MW)
    =#
    @testset "energy_market_can_increase_the_production_level_before_dp_if_ON" begin
        ech = DateTime("2015-01-01T10:00:00")
        next_ech = DateTime("2015-01-01T10:30:00")
        TS = [DateTime("2015-01-01T11:00:00")]
        context, market, firmness = create_instance(ech, next_ech, TS[1],
                                                    22., 0., 33.,
                                                    )
        expected_firmness = PSCOPF.Firmness(
                SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.DECIDED),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.DECIDED), ),
                SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                            "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE), )
                )
        @test firmness == expected_firmness

        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        # prod_1_1 is OFF cause it already was OFF:
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        # prod_1_2 was ON, it can change it's production : 33 => 35.
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
    end


    #=
    ECH1          ECH2            ECH3*          ECH4                  <---FO--->TS
    |             |               |
    10h           10h30           10h40          10h45                          11h
                  <------------------------------------------DP(prod1)----------->
                                  <--------------------------DP(wind1)----------->
                                                   <---------DP(prod2)----------->

    We can no longer change the production level of a unit once DP is reached
    ech = 10h40
    TS = 11h
    DP(prod_1_1) = 30 mins => 10h30 => we are past DP
    DP(prod_1_2) = 15 mins => 10h45

    We suppose that in the preceding step we decided :
    prod_1_1 : 17. (ON)
    prod_1_2 : 20. (ON)

    In terms of cost, it is cheaper to use prod_1_1.
    But, since we are past DP, we can no longer increase the level of prod_1_1.

    wind_1_1 : 20 MW
        remaining Demand : 35MW
    prod_1_1 : 17. (unchanged)
    prod_1_2 : 18. (decreased from 20MW)
    =#
    @testset "energy_market_cannot_change_the_production_level_after_dp" begin
        # We are past DP[prod_1_1]
        ech = DateTime("2015-01-01T10:40:00")
        next_ech = DateTime("2015-01-01T10:45:00")
        TS = [DateTime("2015-01-01T11:00:00")]
        context, market, firmness = create_instance(ech, next_ech, TS[1],
                                                    18., 17., 20.,
                                                    )

        expected_firmness = PSCOPF.Firmness(
                                    SortedDict("prod_1_1" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                ),
                                    SortedDict("wind_1_1" => SortedDict(TS[1] => PSCOPF.TO_DECIDE),
                                                "prod_1_1" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.FREE),
                                                )
            )
        @test firmness == expected_firmness

        # prod_1_1 has a definitive value since we are past the DP

        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        # we are past DP[prod_1_] => we can't change the level of prod_1_1
        # prod_1_1 : 17.
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test 17. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        # prod_1_2 : 20 => 18.
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 18. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
    end

    #=
    at ECH2, we are at wind_1_1's DP,
    but the decision made will still be by scenario cause it depends on the uncertainties
    limitables's DP only applies to limitation decisions (=> applied in TSO)
    ECH1          ECH2            ECH3*          ECH4                  <---FO--->TS
    |             |               |
    10h           10h30           10h40          10h45                          11h
                  <------------------------------------------DP(prod1)----------->
                                  <--------------------------DP(wind1)----------->
                                                   <---------DP(prod2)----------->
    =#
    @testset "energy_market_wind_is_always_uncertain" begin
        # We are past DP[prod_1_1]
        ech = DateTime("2015-01-01T10:40:00")
        next_ech = DateTime("2015-01-01T10:45:00")
        TS = [DateTime("2015-01-01T11:00:00")]
        context, market, firmness = create_instance(ech, next_ech, TS[1],
                                                    18., 17., 20.,
                                                    )

        expected_firmness = PSCOPF.Firmness(
                                    SortedDict("prod_1_1" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                ),
                                    SortedDict("wind_1_1" => SortedDict(TS[1] => PSCOPF.TO_DECIDE),
                                                "prod_1_1" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.FREE),
                                                )
            )
        @test firmness == expected_firmness

        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test !PSCOPF.is_definitive(PSCOPF.get_prod_uncertain_value(context.market_schedule, "wind_1_1", TS[1]))
    end

    #=
    ECH1*         ECH2            ECH3           ECH4                  <---FO--->TS
    |             |               |
    10h           10h30           10h40          10h45                          11h
                  <------------------------------------------DP(prod1)----------->
                                  <--------------------------DP(wind1)----------->
                                                   <---------DP(prod2)----------->
    =#
    @testset "energy_market_cannot_change_the_production_level_before_dp_if_unit_is_off_and_past_DMO" begin
        ech = DateTime("2015-01-01T10:00:00")
        next_ech = DateTime("2015-01-01T10:30:00")
        TS = [DateTime("2015-01-01T11:00:00")]
        context, market, firmness = create_instance(ech, next_ech, TS[1],
                                                    35., 0., 25.,
                                                    )

        expected_firmness = PSCOPF.Firmness(
                                    SortedDict("prod_1_1" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.DECIDED),
                                                ),
                                    SortedDict("wind_1_1" => SortedDict(TS[1] => PSCOPF.FREE),
                                                "prod_1_1" => SortedDict(TS[1] => PSCOPF.FREE),
                                                "prod_1_2" => SortedDict(TS[1] => PSCOPF.FREE),
                                                )
            )
        @test firmness == expected_firmness

        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        # We are past DMO => unit state was already fixed
        # prod_1_1 : OFF
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        # prod_1_2 : 25 => 35.
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
    end

end
