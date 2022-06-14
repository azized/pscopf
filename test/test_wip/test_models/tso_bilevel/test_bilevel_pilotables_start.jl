using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_pilotables_start" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    function create_instance(load_1, load_2,
                            market_decision_prod1, market_decision_prod2,
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 200.,
                                                1000., 10.,
                                                Dates.Second(10), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                20., 200.,
                                                5000., 50.,
                                                Dates.Second(10), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", load_1)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", load_2)

        gen_initial_state = SortedDict{String,PSCOPF.GeneratorState}(
                                "prod_1_1" => PSCOPF.ON,
                                "prod_2_1" => PSCOPF.ON,
                            )

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    gen_initial_state,
                                    uncertainties, nothing, logs)
        if !ismissing(market_decision_prod1)
            PSCOPF.set_commitment_value!(context.market_schedule, "prod_1_1", TS[1], "S1", market_decision_prod1)
        end
        if !ismissing(market_decision_prod2)
            PSCOPF.set_commitment_value!(context.market_schedule, "prod_2_1", TS[1], "S1", market_decision_prod2)
        end

        return context
    end

    @testset "importance_of_input_unit_states_for_tso_bilevel" begin
        #=
        These tests show the importance of the preceding ech decisions as it will affect the TSO decisions
        Here, we illustrate different starting cases.
        Since the preceding modules may have been launched on different uncertainties,
         we can consider arbitrary starting scenarios.
        This proves the importance of the hypothesis that a TSOBilevel model considers a balanced state as input.

        The current model does not allow the TSO to make efficient starting decisions.
         In addition, its limitation decisions affect the performance of the markets.

        The only thing that changes in the testcases is the preceding decision on the unit's start state.
        (these are the decisions at ech-1 for ts and not the generator initial state at ts-1)
        =#


        #=
        TS: [11h]
        S: [S1]
                            bus 1                   bus 2
                            |                      |
        (limitable) prod_1_1|       "1_2"          |prod_2_1
        Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
        Csta=1k, Cprop=10   |                      |Csta=5k, Cprop=50
                            |----------------------|
                            |         35           |
                            |                      |
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 30              |                      | S1: 30

        Suppose market started prod_1_1 and prod_2_1 at the preceding ech.

        There are no RSO constraints
        => TSO does not have to undertake an action
        => TSO keeps the started units at their bounds
        => prod_1_1 : 20 -> 200
           prod_2_1 : 20 -> 200

        Bilevel Market decision would be :
            prod_1_1 : 40 (economically preferred)
            prod_2_1 : 20 (to respect imposed minimum)

        However, the next market module decision without the limitation constraints would have been:
            prod_1_1 : 60 (economically preferred)
            prod_2_1 : 0 (shut down)
            => the limitation constraints lead to higher costs for no RSO reason

        FIXME? we may want the TSO to shutdown one of the units
        =#
        @testset "tso_avoids_shutting_down_units" begin

            context = create_instance(30., 30.,
                                    PSCOPF.ON, PSCOPF.ON,
                                    35.,"start_imp_test")

            tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #Market EOD constraints are OK
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #TSO does not need to bound pilotables because there is no risk of breaking RSO constraints
            #prod_1_1 is not bound
            @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
            #prod_2_1 is not bound
            @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

            #Market chooses the levels respecting the bounds for pilotable production
            @test 40. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
            @test 20. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

            @test - 1e-09 < objective_value(result.upper.model) < 1e-09

        end

        #=
        TS: [11h]
        S: [S1]
                            bus 1                   bus 2
                            |                      |
        (limitable) prod_1_1|                      |prod_2_1
        Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
        Csta=1k, Cprop=10   |        "1_2"         |Csta=5k, Cprop=50
                            |----------------------|
                            |         35           |
                            |                      |
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 30              |                      | S1: 30

        Suppose market did not start any of prod_1_1 nor prod_2_1 (first TSO launch ?)
        The EOD constraint would be violated
        but there are no RSO constraints
        => TSO does not have to undertake any actions cause EOD is not his responsibility
        => TSO keeps the units shut down

        COST for TSO : 0. (no limitations cause units were already shutdown, no LoL for RSO reasons)
        However, if the TSO had started any of the units he would have fixed its minimum production,
        Thus, he would have paid for it.

        PROBLEM : the following market will not be able to start the units due to the impositions
        => prod_1_1 : 0
           prod_2_1 : 0
           But the market will need to cut all the conso
           market cost : LoL(60MW)

        However, the next market module decision without the limitation constraints would have been:
            prod_1_1 : 60 (economically preferred)
            prod_2_1 : 0 (shut down)
            => the limitation constraints lead to LoL for no reason

        FIXME?
        we want TSO to start prod_1_1 without limiations (i.e 20-200) => costs 20. for tso
        =#
        @testset "tso_bilevel_does_not_start_units_for_eod_reasons_from_off" begin
            #no_risk_of_breaking_rso_constraint_off_off

            context = create_instance(30., 30.,
                                    PSCOPF.OFF, PSCOPF.OFF,
                                    35.,"start_imp_test")

            tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            # Solution status
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #Market
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #TSO impositions
            #prod_1_1
            @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
            #prod_2_1
            @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
            @test value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"]) < 1e-09

            #Market chooses the levels respecting the bounds for pilotable production
            @test 60. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
            @test value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"]) < 1e-09

            @test value(PSCOPF.get_upper_obj_expr(result)) ≈ 20. * 10. #imposition cost of prod_1_1
            @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (60 * 10.)

        end

        #=
        same as above :

        TS: [11h]
        S: [S1]
                            bus 1                   bus 2
                            |                      |
        (limitable) prod_1_1|       "1_2"          |prod_2_1
        Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
        Csta=1k, Cprop=10   |                      |Csta=5k, Cprop=50
                            |----------------------|
                            |         35           |
                            |                      |
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 30              |                      | S1: 30

        Suppose we still haven't decided the state of prod_1_1 nor prod_2_1 (first TSO launch ?)
        This is similar to the preceding since we have the same cost expressions for missing or OFF units

        The EOD constraint would be violated
        but there are no RSO constraints
        => TSO does not have to undertake an action cause EOD is not his responsibility
        => he keeps the units shut down
        PROBLEM : the following market will not be able to start the units due to the impositions
        => prod_1_1 : 0
           prod_2_1 : 0
           These cost 0.
           But the market will need to cut all the conso

        If the TSO had started any of the units he would have fixed its minimum production,
        Thus, he would have paid for it.

        FIXME? issue if first tso bilevel (mode 3)
        =#
        @testset "tso_bilevel_does_not_start_units_for_eod_reasons_from_missing" begin
            # no_risk_of_breaking_rso_constraint_missing_missing
            context = create_instance(30., 30.,
                                    missing, missing,
                                    35.,"start_imp_test")

            tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            # Solution status
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #Market
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #TSO impositions
            #prod_1_1
            @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
            #prod_2_1
            @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
            @test value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"]) < 1e-09

            #Market chooses the levels respecting the bounds for pilotable production
            @test 60. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
            @test value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"]) < 1e-09

            @test value(PSCOPF.get_upper_obj_expr(result)) ≈ 20. *10. #imposition of prod_1_1
            @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (60 * 10.)

        end

        #=
        TS: [11h]
        S: [S1]
                            bus 1                   bus 2
                            |                      |
        (limitable) prod_1_1|       "1_2"          |prod_2_1
        Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
        Csta=1k, Cprop=10   |                      |Csta=5k, Cprop=50
                            |----------------------|
                            |         35           |
                            |                      |
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 30              |                      | S1: 30

        Suppose market only started prod_1_1
        There are no RSO constraints
        => TSO does not have to undertake an action
        => he keeps the unit prod_2_1 shutdown
        => prod_1_1 : 20 -> 200 : cost 0.
           prod_2_1 : 0 : cost 0.
           TSO cost : 0.

        The Bilevel market will use prod_1_1 to produce 60MW

        The same behaviour will be observed by the following market module
        =#
        @testset "no_risk_of_breaking_rso_constraint_on_off" begin

            context = create_instance(30., 30.,
                                    PSCOPF.ON, PSCOPF.OFF,
                                    35.,"start_imp_test")

            tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #Market EOD constraints are OK
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #TSO impositions
            #prod_1_1
            @test 20. ≈value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
            #prod_2_1
            @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
            @test value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"]) < 1e-09

            #Market chooses the levels respecting the bounds for pilotable production
            @test 60. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
            @test value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"]) < 1e-09

            @test - 1e-09 < objective_value(result.upper.model) < 1e-09

        end

        #=
        TS: [11h]
        S: [S1]
                            bus 1                   bus 2
                            |                      |
        (limitable) prod_1_1|       "1_2"          |prod_2_1
        Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
        Csta=1k, Cprop=10   |                      |Csta=5k, Cprop=50
                            |----------------------|
                            |         35           |
                            |                      |
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 30              |                      | S1: 30

        Suppose market only started prod_2_1
        There are no RSO constraints
        => TSO does not have to undertake an action
        => TSO keeps the unit prod_1_1 shutdown
        => prod_1_1 : 0 : cost 0.
           prod_2_1 : 20 -> 200 : cost 0.
           These cost 0.

        If TSO had imposed prod_1_1's starting, TSO would have paid : 20
        =#
        @testset "tso_follows_preceding_market_starts" begin

            context = create_instance(30., 30.,
                                    PSCOPF.OFF, PSCOPF.ON,
                                    35.,"start_imp_test")

            tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #Market EOD constraints are OK
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #TSO does not need to bound pilotables
            #prod_1_1 is not bound
            @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
            @test value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]) < 1e-09
            #prod_2_1 is not bound
            @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
            @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

            #Market chooses the levels respecting the bounds for pilotable production
            @test value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"]) < 1e-09
            @test 60. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

            @test objective_value(result.upper.model)  < 1e-09

        end

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 100

    Market chooses prod_1_1 cause cheaper than prod_2_1
    => a flow of 100. => RSO problem

    TSO needs to impose using prod_2_1. two options:

    option 1 : impose prod_1_1 to produce at most 10+35=45 => P_1_1 in [0,45]
        If the TSO chooses to shut down prod_1_1 entirely i.e. [O, O] it will pay :
            => cost (unit was ON): (0-20) + (200-0) = -20 + 200 = 180
        If the TSO chooses to allow using prod_1_1 i.e. [20, 45] it will pay :
            => cost (unit was ON): (20-20) + (200-45) = 155
    (*) If the TSO chooses to allow using prod_1_1 freely i.e. [20, 200] and use option 2 it will pay :
            => cost (unit was ON): (20-20) + (200-200) = 0.    + cost(option 2)

    option 2 : impose prod_2_1 to produce at least 100-35=65 => P_2_1 in [65,200] => cost (unit was OFF): 65

    Imposing one of the preceding conditions guarantees RSO constraint
        => TSO will impose the cheaper one => option 1

    FIXME Note: this might not be what we want, the TSO could have shut bus1 completely and only used bus2
    =#
    @testset "rso_constraint_requires_setting_pilotables_bounds_min" begin

        context = create_instance(10., 100.,
                                PSCOPF.ON, PSCOPF.OFF,
                                35.,"pilotable_test")

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test 65. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels sets the bounds for pilotable production
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 65. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test 65. ≈ objective_value(result.upper.model)

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 300

    Suppose market starts prod_1_1 and prod_2_1
        because prod_1_1 is cheaper, it will be used at max capa i.e. 200.
        This will cause an RSO constraint on branch 1_2

    TSO needs to solve the RSO constraint by
    option 1 : limit using prod_1_1 to 10+35=45
        => limitng prod_1_1 at [20, 45]
        prod_1_1 : 20 -> 45, cost : 200-45 = 155
        prod_2_1 : 20 -> 200, cost : 0.
        => costs 155
    option 2 : cut conso on bus 2 (expensive)

    option 1 is adopted.
    This however will cause an EOD disbalance in market
     cause demand=310, capacity=245
     => market needs to loss_of_load by 65
    =#
    @testset "rso_constraint_requires_setting_pilotables_bounds_max" begin

        context = create_instance(10., 300.,
                                PSCOPF.ON, PSCOPF.ON,
                                35.,"pilotable_test")

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 65. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        # RSO locates cut conso on bus 2 to assure flow limit
        @test value(result.upper.lol_model.p_loss_of_load["bus_1", TS[1],"S1"]) < 1e-09
        @test 65. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2", TS[1],"S1"])

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
        @test 45. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels sets the bounds for pilotable production
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test objective_value(result.upper.model) ≈ ( 155.
                                                        + 65. *tso.configs.TSO_LOL_PENALTY #market LoL
                                                    )

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=20, Pmax=200   |                      |Pmin=20, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         500          |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 400           |                      | S1: 20

    No risk of RSO constraint breaking, but EOD cannot be satisfied
        => the market needs to cut conso
        => the TSO will locate the bus where we cut
    =#
    @testset "eod_problem_requires_cutting_conso" begin

        context = create_instance(400., 20.,
                                PSCOPF.ON, PSCOPF.ON,
                                500.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=false))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints can be respected without capping or losing load OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #Market needs to cut conso to assure EOD constraint
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 20. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        # RSO locates cut conso : here any combination will do
        @test 20. ≈ ( value(result.upper.lol_model.p_loss_of_load["bus_1", TS[1],"S1"])
                    + value(result.upper.lol_model.p_loss_of_load["bus_2", TS[1],"S1"]) )

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels sets the bounds for pilotable production
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

    end


end
