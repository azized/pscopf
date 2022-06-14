using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_pilotables_scenarios" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")

    #=
    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |wind_2_1
    Pmin=0, Pmax=200    |                      |Pmin=0, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------| S1: 30
                        |         35           | S2: 70
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 0             |                      | S1: 100
      S2: 0             |                      | S2: 100


      In scenario 1,
          we can use all the limitable power (30MW),
          this leaves a residual demand of 70MW
          prod_1_1 can only supply 35MW due to RSO constraints
          => a Lol of 35MW
      In scenario 2,
          we can use all the limitable power (70MW),
          this leaves a residual demand of 30MW
          prod_1_1 can supply the required 30MW

      The TSO will need to take an action to limit production of prod_1_1 especially in scenario S1.
    =#
    function create_instance(limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., 10.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "wind_2_1", PSCOPF.Networks.LIMITABLE,
                                                0., 200.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_2_1", TS[1], "S1", 30.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_2_1", TS[1], "S2", 70.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", 0.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S2", 0.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S2", 100.)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, nothing, logs)

        return context
    end

    @testset "different_impositions_by_scenario" begin

        context = create_instance(35.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL=false,
                                LINK_SCENARIOS_PILOTABLE_ON=false))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution has slack
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test_broken 35. ≈ value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) #LoL is for RSO reason but is only noticed in market

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"])

        #wind_2_1
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S1"]) > 30. - 1e-09
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S2"]) > 70. - 1e-09

        #Market chooses the levels sets the bounds for pilotable production
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 30. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"])

        #tso distributes limitable power
        @test 30. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S1"])
        @test 70. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S2"])

    end

    @testset "link_scenarios_in_tso" begin

        context = create_instance(35.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL=true,
                                LINK_SCENARIOS_PILOTABLE_ON=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution has slack
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test_broken 35. ≈ value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) #LoL is for RSO reason but is only noticed in market

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1
        @test ( value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
                ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) )
        @test ( value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
                ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"]) )

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) < 1e-09
        @test 35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"])

        #wind_2_1
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S1"]) > 30. - 1e-09
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S2"]) > 70. - 1e-09

        #Market chooses the levels sets the bounds for pilotable production
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 30. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"])

        #tso distributes limitable power
        @test 30. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S1"])
        @test 70. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S2"])
    end

    #=
    This is not what is illustrated here but keep in mind !

    I'm afraid the TSO can cheat the market when
    LINK_SCENARIOS_PILOTABLE_LEVEL = false
    and
    LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET = true :
    The TSO can penalize one scenario (leaving the others free which costs less)
    but all the market scenarios will obey to the most strict scenario
    since scenarios are linked in the market.
    =#
    @testset "link_scenarios_in_market" begin

        context = create_instance(35.,)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution has slack
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test_broken 35. ≈ value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) #LoL is for RSO reason but is only noticed in market

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]) > 35. - 1e-09

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) < 1e-09
        @test value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"]) > 35. - 1e-09

        @test ( (35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]))
                || (35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"])) )
        @test ( (200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]))
                || (200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"])) )

        #wind_2_1
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S1"]) > 30. - 1e-09
        @test 65. ≈ value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S2"])

        #Market chooses the levels sets the bounds for pilotable production
        @test ( value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
                ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"]) )
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"])

        #tso distributes limitable power
        @test 30. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S1"])
        @test 65. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S2"])

    end

    @testset "link_scenarios" begin

        context = create_instance(35.,)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(LINK_SCENARIOS_PILOTABLE_LEVEL_MARKET=true,
                                                        LINK_SCENARIOS_PILOTABLE_LEVEL=true,
                                                        LINK_SCENARIOS_PILOTABLE_ON=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution has slack
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test_broken 35. ≈ value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) #LoL is for RSO reason but is only noticed in market

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1
        @test ( value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
                ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) )
        @test ( value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
                ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"]) )

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])

        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S2"]) < 1e-09
        @test 35. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S2"])

        #wind_2_1
        @test value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S1"]) > 30. - 1e-09
        @test 65. ≈ value(result.upper.limitable_model.p_limit["wind_2_1",TS[1],"S2"])

        #Market chooses the levels sets the bounds for pilotable production
        @test ( value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
                ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"]) )
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 35. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S2"])

        #tso distributes limitable power
        @test 30. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S1"])
        @test 65. ≈ value(result.upper.limitable_model.p_injected["wind_2_1",TS[1],"S2"])

    end


end
