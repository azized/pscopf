using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_limitables" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    function create_instance(load_1, load_2,
                            wind_1,
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "wind_2_1", PSCOPF.Networks.LIMITABLE,
        #                                         0., 100.,
        #                                         0., 1.,
        #                                         Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], "S1", wind_1)
        # PSCOPF.add_uncertainty!(uncertainties, ech, "wind_2_1", TS[1], "S1", wind_2)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", load_1)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", load_2)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, nothing, logs)

        return context
    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|                      |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |        "1_2"         |
      S1: 40            |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 30

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
    =#
    @testset "no_problem" begin

        context = create_instance(10., 30.,
                                40.,
                                35.)

        tso = PSCOPF.TSOBilevel()
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

        #no limitation
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09

        #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) < 1e-09
        @test value(PSCOPF.get_lower_obj_expr(result)) < 1e-09

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|                      |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |         "1_2"        |
      S1: 100           |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 30

    Available limitable production : 100
    Demand : 40

    Market :
        The market needs to satisfy EOD constraint => only option is to cap 60MW of limitable production.
            e = p_capping = 60.
            lol = p_loss_of_load = 0.

    The TSO locates the capping to avoid RSO constraints (even if there is no risk of violating RSO here)
        p_capping[wind_1_1]=60, plim[wind1]=40, p_injected[wind_1_1]=40

    The TSO takes a limitation action:
      e_min = p_global_capping = 0.
      lol_min = p_global_loss_of_load = 0.
      plim[wind1] = 40 (due to locating capping)

    TSO cost : 0.001 (one limitation)
    Market cost : 60 (capping)

    NOTE: here, there is only one limitable,
     normally the TSO should not have to limit the production of wind_1_1
     cause the market does not have other options
    FIXME? : TSO limits for EOD reasons.
    =#
    @testset "eod_overproduction_requires_market_capping_and_tso_limitation" begin
        context = create_instance(10., 30.,
                                100.,
                                35.,
                                "test_li")

        tso = PSCOPF.TSOBilevel()
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        #TSO :
        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) <= 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) <= 1e-09

        # TSO needs to limit because it is him who locates capping
        @test 40. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"])
        @test 60. ≈ ( value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"]))

        #Market :
        #Market caps ENR for EOD reasons
        @test 60. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S1"])
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) <= 1e-09

        #costs
        @test 0.001 == (1*tso.configs.TSO_LIMIT_PENALTY) #one limitation
        @test ≈(0.001, value(PSCOPF.get_upper_obj_expr(result)), atol=1e-04) #due to the 1e5 objective coeff on global LoL
        @test (60. *tso.configs.MARKET_CAPPING_COST) ≈ value(PSCOPF.get_lower_obj_expr(result))

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 40            |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 100           |                      | S1: 30

    Available limitable production : 40
    Demand : 130

    Market :
        The market needs to satisfy EOD => reduces the consumption to 40
            e = p_capping = 0.
            lol = p_loss_of_load = 90.

    The TSO locates the LoL to avoid RSO constraints
        Here, there is no risk of violating the RSO constraint => any option will do:
            p_loss_of_load[bus_1]+p_loss_of_load[bus_2] = 90

    The TSO takes no further actions
      e_min = p_global_capping = 0.
      lol_min = p_global_loss_of_load = 0.
      no limitation

    TSO cost : 0
    Market cost : 90*1e5 (capping)

    FIXME? : cost LoL locations differently
    =#
    @testset "eod_overload_requires_market_lol" begin
        context = create_instance(100., 30.,
                                40.,
                                35.)

        tso = PSCOPF.TSOBilevel()
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Market needed to cut conso
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO :
        #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) <= 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) <= 1e-09

        #no limitation
        @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
        @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"]) < 1e-09

        #Market :
        #Market cuts conso for EOD reasons
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 90. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        #TSO distributes the cut conso (arbitrarily) while assuring RSO constraint
        @test 90. ≈ ( value(result.upper.lol_model.p_loss_of_load["bus_1", TS[1],"S1"])
                    + value(result.upper.lol_model.p_loss_of_load["bus_2", TS[1],"S1"]) )

        #costs
        @test (90*tso.configs.TSO_LOL_PENALTY) ≈ value(PSCOPF.get_upper_obj_expr(result))
        @test (90*tso.configs.MARKET_LOL_PENALTY) ≈ value(PSCOPF.get_lower_obj_expr(result))

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 50            |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 40


    Available limitable production : 50
    Demand : 50

    Market :
        The EOD is satisfied, the market would simply use the limitable production to satisfy the demand.
        => this will cause a flow[branch_1_2] = 40 > 35  ==> violates the RSO constraint

    The TSO needs to impose some additional constraints to oblige the market to respect the RSO constraint:
        option 1 : impose a capping of 5MW, limit wind_1_1's production to 45MW
        => cost : 5 + 0.001
            This will oblige the market to cut conso by 5.
            Since the TSO locates the cut conso, he can assure cutting conso on bus2 to assure the RSO constraint
        option 2 : impose reducing the consumption of bus2 by 5MW, limit wind_1_1's production to 45MW
            bus 2 will only require 35MW which assures RSO constraints
            TSO will need to limit the prod of wind_1_1
        => cost : 5e05 + 0.001

    TSO adopts the cheaper option : option 1
        e_min = p_global_capping = 5.
        lol_min = p_global_loss_of_load = 0.
        plimit[wind_1_1] = 45

    Market
    The market is now constrained by the TSO decisions:
        e = p_capping = 5. (due to TSO decisions)
        lol = p_loss_of_load = 5. (due to EOD)

    The TSO locates the capping and LoL:
        p_capping[wind_1_1]=5
        p_loss_of_load[bus_1]=0, p_loss_of_load[bus_2]=5

    TSO cost : 5.001
    Market cost : 5 + 5*1e5

    #FIXME? TSO obliges cuts for RSO reasons
    #normally both capping and lol are due to RSO constraints
    #however the TSO simply tells the market to cap some MW which only costs 5€
    =#
    @testset "rso_problem" begin
        context = create_instance(10., 40.,
                                50.,
                                35.)

        tso = PSCOPF.TSOBilevel()
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Market needed to cut conso
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO :
        #TSO imposes caping 5MW
        @test 5. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S1"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #limit wind_1_1
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"])
        @test 5. ≈ value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"])

        #Market :
        #Market only cuts as required since EOD is OK
        @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S1"])
        @test 5. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        #TSO decides the distribution :
        @test 5. ≈ ( value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"]))
        @test 5. ≈ ( value(result.upper.lol_model.p_loss_of_load["bus_1",TS[1],"S1"])
                    + value(result.upper.lol_model.p_loss_of_load["bus_2",TS[1],"S1"]) )
        @test value(result.upper.lol_model.p_loss_of_load["bus_1",TS[1],"S1"]) < 1e-09
        @test 5. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2",TS[1],"S1"])

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈  (5*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY
                                                            + 5*tso.configs.TSO_LOL_PENALTY )
        @test (5*tso.configs.MARKET_CAPPING_COST + 5*tso.configs.MARKET_LOL_PENALTY) ≈ value(PSCOPF.get_lower_obj_expr(result))
    end

    #=
    same as preceding : c.f. rso_problem

    option 1 : impose a capping of 5MW, limit wind_1_1's production to 45MW
    => default parameters cost : 5 + 0.001
    => cost : 5e03 + 0.001
    option 2 : impose reducing the consumption of bus2 by 5MW, limit wind_1_1's production to 45MW
    => default parameters cost : 5e05 + 0.001
    => cost : 5*3 + 0.001

    with new config, option 2 is more appealing
    =#
    @testset "rso_problem_with_modifed_costs" begin
        context = create_instance(10., 40.,
                                50.,
                                35.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(
                                            TSO_LOL_PENALTY = 3.,
                                            TSO_CAPPING_COST = 1e3
                                ))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Market needed to cut conso
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        #TSO :
        #TSO imposes reducing the consumption
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test 5. ≈ value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"])

        #limit wind_1_1
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"])
        @test 5. ≈ value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"])

        #Market :
        #Market only cuts as required since EOD is OK
        @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S1"])
        @test 5. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        #TSO decides the distribution :
        @test 5. ≈ ( value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"]))
        @test 5. ≈ ( value(result.upper.lol_model.p_loss_of_load["bus_1",TS[1],"S1"])
                    + value(result.upper.lol_model.p_loss_of_load["bus_2",TS[1],"S1"]) )
        @test value(result.upper.lol_model.p_loss_of_load["bus_1",TS[1],"S1"]) < 1e-09
        @test 5. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2",TS[1],"S1"])

        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ 30.001 ≈ (5*tso.configs.TSO_LOL_PENALTY + #due to TSO LoL
                                                                    1*tso.configs.TSO_LIMIT_PENALTY +
                                                                    5*tso.configs.TSO_LOL_PENALTY) #due to market LoL
        @test (5+5e05) == (5*tso.configs.MARKET_CAPPING_COST + 5*tso.configs.MARKET_LOL_PENALTY) ≈ value(PSCOPF.get_lower_obj_expr(result))
    end

end
