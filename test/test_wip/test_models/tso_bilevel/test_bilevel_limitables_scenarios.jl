using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_limitables_scenarios" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    scenarios = ["S1", "S2"]
    function create_instance(scenarios,
                            load_1::Vector{Float64}, load_2::Vector{Float64}, wind_1::Vector{Float64},
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        for (s_index, s) in enumerate(scenarios)
            PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], s, load_1[s_index])
            PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], s, load_2[s_index])
            PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], s, wind_1[s_index])
        end

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, nothing, logs)

        return context
    end

    #=
    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 40            |----------------------|
      S2: 60            |         35           |
                        |                      |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 30
      S2: 25            |                      | S2: 35

    S1 :
        Available limitable production : 40
        Demand : 40

        Market : use wind_1_1 at 40.
        This satisfies the EOD constraint
            e = p_capping = 0.
            lol = p_loss_of_load = 0.
        This does not induce any RSO constraint violation

        The TSO does not need to take any actions :
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            no limitation

        TSO cost : 0
        Market cost : 0

    S2 :
        Available limitable production : 60
        Demand : 60

        Market : use wind_1_1 at 60.
        This satisfies the EOD constraint
            e = p_capping = 0.
            lol = p_loss_of_load = 0.
        This does not induce any RSO constraint violation

        The TSO does not need to take any actions :
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            no limitation

        TSO cost : 0
        Market cost : 0
    =#
    @testset "no_problem" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [10., 25.], #load1
                                    [30., 35.], #load2
                                    [40., 60.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        tso.configs.LINK_SCENARIOS_LIMIT = false

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

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        # no limitations : Limitations are by scenario : p_limit can be different for scenario 1 and 2
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"]) > 60. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"]) < 1e-09

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) < 1e-09
        @test value(PSCOPF.get_lower_obj_expr(result)) < 1e-09

    end

    #=
    previous + link limits : no change cause there was no limitations.
    =#
    @testset "no_problem_link_scenarios_limit" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [10., 25.], #load1
                                    [30., 35.], #load2
                                    [40., 60.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        #default, tso.configs.LINK_SCENARIOS_LIMIT = true

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

        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S2"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        # Limitations are by scenario
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"]) > 60. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"]) < 1e-09
        @test !PSCOPF.is_different( value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]),
                                value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"]) )

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) < 1e-09
        @test value(PSCOPF.get_lower_obj_expr(result)) < 1e-09

    end

    #=
    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 40            |----------------------|
      S2: 60            |         35           |
                        |                      |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 30
      S2: 10            |                      | S2: 50

    S1 :
        Available limitable production : 40
        Demand : 40

        Market : use wind_1_1 at 40.
        This satisfies the EOD constraint
            e = p_capping = 0.
            lol = p_loss_of_load = 0.
        This does not induce any RSO constraint violation

        The TSO does not need to take any actions :
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            no limitation

        TSO cost : 0
        Market cost : 0

    S2 :
        Available limitable production : 60
        Demand : 60

        Market :
            EOD can be satisfied. The market would use wind_1_1 at 60. But, this will cause a RSO constraint.
            => TSO needs to take an action to prevent this.

        TSO has two options :
            option 1 : impose a capping of 15MW, limit wind_1_1's production to 45MW
            option 2 : impose reducing the consumption of bus2 by 15MW, limit wind_1_1's production to 45MW
        He chooses the cheaper : option 1
            e_min = p_global_capping = 15.
            lol_min = p_global_loss_of_load = 0.
            plimit[wind_1_1] = 45

        Market
        The market is now constrained by the TSO decisions:
            e = p_capping = 15. (due to TSO decisions)
            lol = p_loss_of_load = 15. (due to EOD)

        The TSO locates the capping and LoL:
            p_capping[wind_1_1]=15
            p_loss_of_load[bus_1]=0, p_loss_of_load[bus_2]=15

        TSO cost : 15.001
        Market cost : 15 + 15*1e5


    Total : S1 + S2
        TSO cost : 15.001
        Market cost : 15 + 15*1e5

    =#
    @testset "eod_problem_needs_capping_in_s2" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [10., 10.], #load1
                                    [30., 50.], #load2
                                    [40., 60.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        tso.configs.LINK_SCENARIOS_LIMIT = false

        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK in S1, but not in S2
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 15. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK, but there would be a disbalance in S2 due to RSO action
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 15. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test 15. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"])

        # need to limit in S2
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"])
        @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"])

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (
                    0. #S1
                    + 15*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 15*tso.configs.TSO_LOL_PENALTY#S2
                    )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (
                    0. #S1
                    + 15*tso.configs.MARKET_CAPPING_COST + 15*tso.configs.MARKET_LOL_PENALTY #S1
                    )

    end

    #=
    c.f. eod_problem_needs_capping_in_s2

    previous + link limits : no changes because the limit of S2 (limited scenario)
    can also be considered as a limit for S1 since it's higher than the uncertainty.
    =#
    @testset "eod_problem_needs_capping_in_s2_link_scenarios" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [10., 10.], #load1
                                    [30., 50.], #load2
                                    [40., 60.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        #Default, tso.configs.LINK_SCENARIOS_LIMIT = true

        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO needs to cap prod at S2 due to RSO constraints
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 15. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are disbalanced due to TSO limiting
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 15. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test 15. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"])

        # Limitations are by scenario
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"])
        @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"])
        @test !PSCOPF.is_different( value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]),
                                value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"]) )

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (
                    0. #S1 (no limitation)
                    + 15*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 15*tso.configs.TSO_LOL_PENALTY #S2
                    )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (
                    0. #S1
                    + 15*tso.configs.MARKET_CAPPING_COST + 15*tso.configs.MARKET_LOL_PENALTY #S1
                    )

    end

    #=
    TS: [11h]
    S: [S1, S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 60            |----------------------|
      S2: 65            |         35           |
                        |                      |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 50            |                      | S1: 10
      S2: 10            |                      | S2: 55

    S1 :
        Available limitable production : 60
        Demand : 60

        Market : use wind_1_1 at 60.
        This satisfies the EOD constraint
            e = p_capping = 0.
            lol = p_loss_of_load = 0.
        This does not induce any RSO constraint violation

        The TSO does not need to take any actions :
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            no limitation

        TSO cost : 0
        Market cost : 0

    S2 :
        Available limitable production : 65
        Demand : 65

        Market :
            EOD can be satisfied. The market would use wind_1_1 at 65.
            But, this will cause a RSO constraint : a flow of 55 (>35) on branch 1->2
            => TSO needs to take an action to prevent this.

        TSO has two options :
            option 1 : impose a capping of 20MW, limit wind_1_1's production to 45MW
            option 2 : impose reducing the consumption of bus2 by 20MW, limit wind_1_1's production to 45MW
        He chooses the cheaper : option 1
            e_min = p_global_capping = 20.
            lol_min = p_global_loss_of_load = 0.
            plimit[wind_1_1] = 45

        Market
        The market is now constrained by the TSO decisions:
            e = p_capping = 20. (due to TSO decisions)
            lol = p_loss_of_load = 20. (due to EOD)

        The TSO locates the capping and LoL:
            p_capping[wind_1_1]=20
            p_loss_of_load[bus_1]=0, p_loss_of_load[bus_2]=20

        TSO cost : 20.001
        Market cost : 20 + 20*1e5

    The two scenarios decisions are independant : OK
    =#
    @testset "limiting_a_scenario_does_not_affect_the_other" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [50., 10.], #load1
                                    [10., 55.], #load2
                                    [60., 65.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        tso.configs.LINK_SCENARIOS_LIMIT = false

        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO RSO constraints are OK in S1, but not in S2
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 20. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are OK, but there would be a disbalance in S2 due to RSO action
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 20. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test 20. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"])

        # need to limit in S2
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 60. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"])
        @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"])

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (
                    0. #S1 (no limitation)
                    + 20*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 20*tso.configs.TSO_LOL_PENALTY #S2
                    )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (
                    0. #S1
                    + 20*tso.configs.MARKET_CAPPING_COST + 20*tso.configs.MARKET_LOL_PENALTY #S1
                    )

    end

    #=
    previous + link limits

    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 60            |----------------------|
      S2: 65            |         35           |
                        |                      |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 50            |                      | S1: 10
      S2: 10            |                      | S2: 55

    This is the same previous instance. However, limitation decisions are common among scenarios.

    Initially, No problem in S1
    In S2, RSO constraint is violated, need to limit flow 1->2 to 35
        => limit wind1 to 45=10+35 (option 1)
        or cut conso on bus2 to 35 (option 2)

    capping prod and limiting wind1 to 45 in S2,
        will oblige to limit wind1 to 45 in S1 as well
        =>  for tso : 15MW capped in S1, 20MW capped in S2
                    no need to impose cutting consumption (lol) cause capping + EOD will oblige the market to do it
                    => cost : 35
            for market : 15MW capped in S1, 20MW capped in S2
                    cut 15MW consumption in S1 and 20 MW consumption in S2
                    => cost 35 + 35*1e4
    =#
    @testset "limiting_a_scenario_affects_the_other_link_scenarios" begin

        context = create_instance(["S1", "S2"], #scenarios
                                    [50., 10.], #load1
                                    [10., 55.], #load2
                                    [60., 65.], #wind_1_1
                                35.)

        tso = PSCOPF.TSOBilevel()
        #Default, tso.configs.LINK_SCENARIOS_LIMIT = true

        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        # TSO needs to cap prod at S2 due to RSO constraints,
        #but limiting is the same across scenarios => capping will be needed to solve S1's limiting as well
        @test 15. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S1"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        @test 20. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S2"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S2"]) < 1e-09

        #Market EOD constraints are disbalanced due to TSO limiting
        @test 15. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S1"])
        @test 15. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        @test 20. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S2"])
        @test 20. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S2"])

        # Limitations are common to both scenarios
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"])
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"])
        @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S2"])
        @test !PSCOPF.is_different( value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]),
                                value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S2"]) )

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (
                    15*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 15*tso.configs.TSO_LOL_PENALTY #S1
                    + 20*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 20*tso.configs.TSO_LOL_PENALTY #S2
                    )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (
                    15*tso.configs.MARKET_CAPPING_COST + 15*tso.configs.MARKET_LOL_PENALTY #S1
                    + 20*tso.configs.MARKET_CAPPING_COST + 20*tso.configs.MARKET_LOL_PENALTY #S2
                    )

    end

end
