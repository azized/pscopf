using PSCOPF

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_tso_start_cost" begin


    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1
                        |
    (pilotable) prod_1_1|          load_1
    Pmin=10, Pmax=100   |S1: 30     S1: 50
    Csta=450, Cprop=10  |S2: 25     S2: 35
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=800, Cprop=15 |
    =#

    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T07:00:00")
    network = PSCOPF.Networks.Network()
    # Buses
    PSCOPF.Networks.add_new_bus!(network, "bus_1")
    # Pilotables
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            450., 10.,
                                            Dates.Second(0), Dates.Second(0))
    prod_1_1 = PSCOPF.safeget_generator(network, "prod_1_1")
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_2", PSCOPF.Networks.PILOTABLE,
                                            10., 100.,
                                            800., 15.,
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
    Market has not set these units' levels
    We only use prod_1_1 cause its prop cost is lower
    =#
    @testset "tso_no_starting_cost_if_units_initially_started" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.ON,
            "prod_1_2" => PSCOPF.ON
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) < 1e-09

        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")

        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    prod_1_1 is OFF, prod_1_2 is ON
    We only use prod_1_2 cause its already started => no start cost => cheaper
    =#
    @testset "tso_no_starting_cost_for_unit_prod_1_2" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.ON
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) < 1e-09
        #high start cost for prod_1_1 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #no start cost for prod_1_2 => used eventhough it's unitary prop_cost is higher
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    prod_1_1 is OFF, prod_1_2 is OFF
    We only use prod_1_1 cause it's cheaper (mainly in terms of start cost)
    =#
    @testset "tso_starting_cost" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 2*PSCOPF.get_start_cost(prod_1_1) ≈ value(result.objective_model.start_cost) # started in 2 scenarios

        #start cost for prod_1_1 prefered to prod_1_2
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    Init : prod_1_1 is OFF, prod_1_2 is OFF
    Market : uses prod_1_2 ; started (ON) and set the prod levels satisfying the EOD
    We use prod_1_2 to follow the market, even if prod_1_1 is cheaper
    =#
    @testset "tso_starting_cost_due_to_following_market" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>30.,"S2"=>25.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>50.,"S2"=>35.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 2*PSCOPF.get_start_cost(prod_1_2) ≈ value(result.objective_model.start_cost) # started in 2 scenarios
        @test value(result.objective_model.deltas) < 1e-09 # followed the market

        #start cost for prod_1_1 prefered to prod_1_2 but will cause high deviation from the market
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 but we follow the market
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    Init : prod_1_1 is OFF, prod_1_2 is OFF
    Market : uses prod_1_2 ; started (ON) but didn't set the prod levels
    => using prod_1_1 or prod_1_2 will have the same effect on deltas
    We use prod_1_1 cause cheaper
    =#
    @testset "tso_starting_cost_deviate_from_market" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test 2*PSCOPF.get_start_cost(prod_1_1) ≈ value(result.objective_model.start_cost) # started in 2 scenarios
        @test (30 + 25 + 50 + 35) ≈ value(result.objective_model.deltas)

        #start cost for prod_1_1 prefered to prod_1_2
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    Init : prod_1_1 is OFF, prod_1_2 is OFF
    Market : definitively starts prod_1_2
    => if we use prod_1_2 we do not pay for starting
    =#
    @testset "tso_no_starting_cost_for_definitive_market_starts" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) <= 1e-09 # starting cost was paid by market

        #start cost for prod_1_1 prefered to prod_1_2
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    Init : prod_1_1 is OFF, prod_1_2 is OFF
    TSO already definitively started prod_1_2
    => if we use prod_1_2 we do not pay for starting
    =#
    @testset "tso_no_starting_cost_for_past_definitive_tso_starts" begin
        # initial generators state
        generators_init_state = SortedDict(
            "prod_1_1" => PSCOPF.OFF,
            "prod_1_2" => PSCOPF.OFF
        )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.OFF,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.OFF,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )
        context.tso_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S2"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )
        firmness_l = PSCOPF.Firmness(
                SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.DECIDED,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.DECIDED) ),
                SortedDict( "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE),
                            "prod_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE,
                                                    Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.FREE), )
                )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness_l,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness_l, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) <= 1e-09 # starting cost was paid by market

        #start cost for prod_1_1 prefered to prod_1_2
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1
                        |
    (pilotable) prod_1_1|          load_1
    Pmin=10, Pmax=100   |S1: 30     S1: 50
    Csta=450, Cprop=10  |S2: 25     S2: 165
                        |
    (pilotable) prod_1_2|
     Pmin=10, Pmax=100  |
     Csta=800, Cprop=15 |

    We have high demand for TS2, S2
    =#

    #=
    both units are ON
    In TS2,S2, we need both units
    but in TS1 one unit is enough
    => if we turn off a unit in TS1, we will need to restart it for TS2
    => it's cheaper to use both units at TS1 non-optimally to get the global optimum
    =#
    @testset "tso_both_units_used_in_S2_to_avoid_restart_cost" begin
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

        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test value(result.objective_model.start_cost) <= 1e-09 # No start cost, in S2, prod_1_2 is used in TS1 to avoid restart

        #start cost for prod_1_1 prefered to prod_1_2
        @test 30. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test 15. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test 50. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test 100. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        #higher start cost for prod_1_2 => not used
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test 10. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test 65. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

    #=
    both units are OFF
    in TS1 one unit is enough => start prod_1_1
    In TS2,S2, we need both units => start prod_1_2
    =#
    @testset "tso_both_units_are_started_when_needed" begin
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
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                        SortedDict("S2"=>PSCOPF.OFF, "S2"=>PSCOPF.OFF))
                            ),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                            Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test (2*450 + PSCOPF.get_start_cost(prod_1_2)) ≈ value(result.objective_model.start_cost)  # prod1 started in 2 scenarios, prod2 in S2

        # S1 : prod_1_1 is cheaper to start => used
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test 30. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[2], "S1") < 1e-09

        # S2 : prod_1_2 is needed to satisfy demand it also has high proportional cost
        # => prod_1_2 is only started when needed
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
        @test 25. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[1], "S2") < 1e-09
        @test 65. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end


    #=
    definitive starting decisions made in
     the tso_schedule and the reference market schedule are gratis for the TSO
    Here, in tso_schedule, we :
        - consider starting prod_1_1  (non-definitive decision) => it is not gratis for market
        - start prod_1_2 (definitive decision) => prod_1_2 has no start cost for market
    =#
    @testset "tso_test_gratis_start" begin
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

        context.market_schedule = PSCOPF.Schedule(PSCOPF.TSO(), Dates.DateTime("2015-01-01T06:00:00"), SortedDict(
                                        "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(missing,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                                            ),
                                        "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                                                    SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.ON))
                                                        ),
                                            SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                        Dates.DateTime("2015-01-01T11:15:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>0.)),
                                                                                        )
                                            )
                                        )
                                )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        @test (450) ≈ value(result.objective_model.start_cost)  # prod1 started in 2 scenarios, prod2 in S2
        @test( value(result.objective_model.prop_cost) ≈
              (   (0.   *10. + 30. * 15) #TS1, S1
                + (0.   *10. + 50. * 15) #TS2, S1
                + (15.  *10. + 10.  * 15) #TS1, S2
                + (100. *10. + 65. * 15) #TS2, S2
              )
        )

        #S1 : need one generator, prod_1_2 can be started gratis => no need to start prod_1_1
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test PSCOPF.OFF == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S1")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S1") < 1e-09
        @test 30. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[1], "S1")
        @test 50. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[2], "S1")
        #S2 : need both generators at TS2, prod_1_2 can be started gratis at TS1
        # => we use prod_1_2 since ts1, to get its free starting (starting it at TS2 would cost less prop_cost but <e would pay for starting)
        # prod_1_1 is use right from TS1 even if not really needed cause we will pay for starting anyways, so at least use it at TS1 to reduce prop_cost
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        @test PSCOPF.ON == PSCOPF.get_commitment_value(context.tso_schedule, "prod_1_2", TS[2], "S2")
        @test 15. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S2")
        @test 10. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[1], "S2")
        @test 100. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[2], "S2")
        @test 65. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_2", TS[2], "S2")

        #reset original value
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:15:00"), "S2", 65.)
    end

end