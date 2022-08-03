using PSCOPF

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_energy_market_start_cost" begin

    LOL_COST = 1e5

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1
                        |
    (pilotable) prod_1_1|          load_1
    Pmin=10, Pmax=100   |S1: 30     S1: 50
    Csta=4.5k, Cprop=10 |S2: 25     S2: 35
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=8k, Cprop=15  |
    =#

    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T07:00:00")
    network = PSCOPF.Networks.Network()
    # Buses
    PSCOPF.Networks.add_new_bus!(network, "bus_1")
    # Pilotables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            4500., 10.,
                                            Dates.Second(0), Dates.Second(0))
    prod_1_1 = PSCOPF.safeget_generator(network, "prod_1_1")
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_2", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            8000., 15.,
                                            Dates.Second(0), Dates.Second(0))
    prod_1_2 = PSCOPF.safeget_generator(network, "prod_1_2")
    # Uncertainties
    uncertainties = PSCOPF.Uncertainties()
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 30.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S1", 50.)
    PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 35.)
    # firmness
    firmness = PSCOPF.Firmness(
                SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE) ),
                SortedDict( "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE), )
                )

    #=
    prod_1_1 is ON, prod_1_2 is ON
    We only use prod_1_1 cause its prop cost is lower
    =#
    @testset "energy_market_no_starting_cost_if_units_initially_started" begin

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.ON,
            "prod_1_2" => PSCOPF.ON
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) < 1e-09

        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")

        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    prod_1_1 is OFF, prod_1_2 is ON
    We only use prod_1_2 cause its already started => no start cost => cheaper
    =#
    @testset "energy_market_no_starting_cost_for_unit_prod_1_2" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.ON
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) < 1e-09
        #high start cost for prod_1_1 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        #no start cost for prod_1_2 => used eventhough it's unitary prop_cost is higher
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    prod_1_1 is OFF, prod_1_2 is OFF
    We only use prod_1_1 cause it's cheaper (mainly in terms of start cost)
    =#
    @testset "energy_market_starting_cost" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 2*PSCOPF.get_start_cost(prod_1_1) ≈ value(result.objective_model.start_cost) # started in 2 scenarios

        @assert PSCOPF.get_start_cost(prod_1_1) < PSCOPF.get_start_cost(prod_1_2)

        #start cost for prod_1_1 prefered to prod_1_2
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1
                        |
    (pilotable) prod_1_1|          load_1
    Pmin=10, Pmax=100   |S1: 30     S1: 50
    Csta=4.5k, Cprop=10 |S2: 25     S2: 165
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=8k, Cprop=15  |

    We have high demand for TS2, S2
    =#

    #=
    both units are ON
    In TS2,S2, we need both units
    but in TS1 one unit is enough
    => if we turn off a unit in TS1, we will need to restart it for TS2
    => it's cheaper to use both units at TS1 non-optimally to get the global optimum
    =#
    @testset "energy_market_both_units_used_in_S2_to_avoid_restart_cost" begin
        #bus2, S2, TS:11h15 : high demand
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 165.)

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.ON,
            "prod_1_2" => PSCOPF.ON
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) <= 1e-09 # No start cost, in S2, prod_1_2 is used in TS1 to avoid restart

        #start cost for prod_1_1 prefered to prod_1_2
        @test 30. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test 15. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 50. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 10. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

    #=
    both units are OFF
    in TS1 one unit is enough => start prod_1_1
    In TS2,S2, we need both units => start prod_1_2
    =#
    @testset "energy_market_both_units_are_started_when_needed" begin
        #bus2, S2, TS:11h15 : high demand
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 165.)

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test (2*PSCOPF.get_start_cost(prod_1_1) + PSCOPF.get_start_cost(prod_1_2)) ≈ value(result.objective_model.start_cost)  # prod1 started in 2 scenarios, prod2 in S2

        # S1 : prod_1_1 is cheaper to start => used
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test 30. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S1") < 1e-09

        # S2 : prod_1_2 is needed to satisfy demand it also has high proportional cost
        # => prod_1_2 is only started when needed
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
        @test 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S2") < 1e-09
        @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end


    #=
    Starting decisions made in the TSO actions are gratis for the EnergyMarket
    Here, in tso_actions, we firmly decided :
        - start prod_1_2 at TS1 (definitive decision) => prod_1_2 has no start cost for market at TS1
    prod_1_1  is not listed in TSOActions => it is not gratis for market (the decision was not firm yet)
    =#
    @testset "energy_market_no_gratis_starts" begin
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST,
                                                                CONSIDER_GRATIS_STARTS = false))
        @assert( PSCOPF.is_tso(market.configs.REF_SCHEDULE_TYPE) )

        #bus2, S2, TS:11h15 : high demand
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 165.)

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)

        firmness_l = PSCOPF.Firmness(
            SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.TO_DECIDE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.TO_DECIDE) ),
            SortedDict( "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE), )
            )

        context.tso_schedule = PSCOPF.Schedule(PSCOPF.TSO(), Dates.DateTime("2015-01-01T06:00:00"), SortedDict(
                                        "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>15.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>100.)),
                                                        ),
                                            ),
                                        "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>10.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>65.)),
                                                        ),
                                            )
                                        )
                                )

        result = PSCOPF.run(market, ech, firmness_l,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness_l, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # S2 requires two units, prod_1_2 must be firm => used in both scenarios (gratis)
        @test (PSCOPF.get_start_cost(prod_1_1) + 2 * PSCOPF.get_start_cost(prod_1_2)) ≈ value(result.objective_model.start_cost) # No cost for starting prod_1_2 at TS1
        @test( value(result.objective_model.prop_cost) ≈
              (   (0.   *10. + 30. * 15) #TS1, S1
                + (0.   *10. + 50. * 15) #TS2, S1
                + (15.  *10. + 10.  * 15) #TS1, S2
                + (100. *10. + 65. * 15) #TS2, S2
              )
        )

        #S1 : need one generator, prod_1_1 is cheaper
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S1") < 1e-09
        @test 30. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        #S2 : need both generators at TS2, prod_1_2 is expensive => only started when needed
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
        @test 15. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test 10. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

    #=
    definitive starting decisions made in
     the reference TSO schedule are not gratis for the EnergyMarket
     (the reference for commitment are the TSOActions)
    Here, in tso_schedule, we :
        - consider starting prod_1_1  (non-definitive decision)
        - start prod_1_2 (definitive decision)
    We still need to pay for both options cause it is the TSO who validates commitments.
    and TSO Actions (here, empty) do not indicate that we started prod_1_2 (nor prod_1_1)
    => market pays for both units
    =#
    @testset "energy_market_does_not_consider_gratis_starts_from_preceding_tso_schedule" begin
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        @assert( PSCOPF.is_tso(market.configs.REF_SCHEDULE_TYPE) )

        #bus2, S2, TS:11h15 : high demand
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 165.)

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)

        firmness_l = PSCOPF.Firmness(
            SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.TO_DECIDE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.TO_DECIDE) ),
            SortedDict( "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE), )
            )

        context.tso_schedule = PSCOPF.Schedule(PSCOPF.TSO(), Dates.DateTime("2015-01-01T06:00:00"), SortedDict(
                                        "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>15.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>100.)),
                                                        ),
                                            ),
                                        "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>10.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>65.)),
                                                        ),
                                            )
                                        )
                                )

        result = PSCOPF.run(market, ech, firmness_l,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness_l, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # S2 requires two units, prod_1_2 must be firm => used in both scenarios
        @test PSCOPF.get_start_cost(prod_1_1) ≈ value(result.objective_model.start_cost)
        @test( value(result.objective_model.prop_cost) ≈
              (   (0.   *10. + 30. * 15) #TS1, S1
                + (0.   *10. + 50. * 15) #TS2, S1
                + (15.  *10. + 10.  * 15) #TS1, S2
                + (100. *10. + 65. * 15) #TS2, S2
              )
        )

        #S1 : need one generator, prod_1_1 is cheaper
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S1") < 1e-09
        @test 30. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        #S2 : need both generators at TS2, prod_1_2 is expensive => only started when needed
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
        @test 15. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test 10. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

    #=
    definitive starting decisions made in
     the reference market schedule are not gratis for the EnergyMarket
     (the reference for commitment are the TSOActions)
    Here, in tso_schedule, we :
        - consider starting prod_1_1  (non-definitive decision)
        - start prod_1_2 (definitive decision)
    We still need to pay for both options cause it is the TSO who validates commitments.
    and TSO Actions (here, empty) do not indicate that we started prod_1_2 (nor prod_1_1)
    => market pays for both units
    =#
    @testset "energy_market_does_not_consider_gratis_starts_from_preceding_market_schedule" begin
        market = PSCOPF.EnergyMarket(PSCOPF.EnergyMarketConfigs(loss_of_load_penalty=LOL_COST))
        @assert( PSCOPF.is_tso(market.configs.REF_SCHEDULE_TYPE) )

        #bus2, S2, TS:11h15 : high demand
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 165.)

        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)

        firmness_l = PSCOPF.Firmness(
            SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.TO_DECIDE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.TO_DECIDE) ),
            SortedDict( "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                        "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE), )
            )

        context.market_schedule = PSCOPF.Schedule(PSCOPF.TSO(), Dates.DateTime("2015-01-01T06:00:00"), SortedDict(
                                        "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>15.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>100.)),
                                                        ),
                                            ),
                                        "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>15.,"S2"=>10.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                                                    SortedDict("S1"=>25.,"S2"=>65.)),
                                                        ),
                                            )
                                        )
                                )

        result = PSCOPF.run(market, ech, firmness_l,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness_l, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # S2 requires two units, prod_1_2 must be firm => used in both scenarios
        @test (PSCOPF.get_start_cost(prod_1_1) + 2 * PSCOPF.get_start_cost(prod_1_2)) ≈ value(result.objective_model.start_cost)
        @test( value(result.objective_model.prop_cost) ≈
              (   (0.   *10. + 30. * 15) #TS1, S1
                + (0.   *10. + 50. * 15) #TS2, S1
                + (15.  *10. + 10.  * 15) #TS1, S2
                + (100. *10. + 65. * 15) #TS2, S2
              )
        )

        #S1 : need one generator, prod_1_1 is cheaper
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S1") < 1e-09
        @test 30. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S1")
        #S2 : need both generators at TS2, prod_1_2 is expensive => only started when needed
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.market_schedule, "prod_1_2", TS[2], "S2")
        @test 15. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[2], "S2")
        @test 10. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[1], "S2")
        @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

end
