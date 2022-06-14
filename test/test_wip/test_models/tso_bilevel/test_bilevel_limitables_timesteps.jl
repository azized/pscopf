using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_limitables_timesteps" begin

    TS = [DateTime("2015-01-01T11:00:00"), DateTime("2015-01-01T11:15:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    function create_instance(timesteps,
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
        for (ts_index, ts) in enumerate(timesteps)
            PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", ts, "S1", load_1[ts_index])
            PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", ts, "S1", load_2[ts_index])
            PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", ts, "S1", wind_1[ts_index])
        end

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, nothing, logs)

        return context
    end


    @testset "balanced_ts1" begin
        #=
        TS: [11h, 11h15]
        S: S1
            TS1       TS2
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
            S1: 40    S1: ? |----------------------|
                            |         35           |
                            |                      |
                            |                      |
                load(bus_1) |                      |load(bus_2)
            S1: 10    S1: ? |                      | S1: 30   S1: ?

        AT 11h :

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

        AT 11h15 : uncertainties values at 11h15 will change in each test.
                These won't affect decisions for 11h.
        =#
        function  test_ts1_values(result, TS)
            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

            #no limitation
            @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"]) > 40. - 1e-09
            @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"]) < 1e-09

            #Market EOD constraints are OK
            @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        end


        #=
        TS: [11h, 11h15]
        S: S1
            TS1       TS2
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
            S1: 40    S1: 60|----------------------|
                            |         35           |
                            |                      |
                            |                      |
                load(bus_1) |                      |load(bus_2)
            S1: 10    S1: 25|                      | S1: 30   S1: 35

        At 11h15:
            Available limitable production : 60
            Demand : 60

            Market : use wind_1_1 at 60.
            This satisfies the EOD constraint
                e = p_capping = 0.
                lol = p_global_loss_of_load = 0.
            This does not induce any RSO constraint violation

            The TSO does not need to take any actions :
                e_min = p_global_capping = 0.
                lol_min = p_global_loss_of_load = 0.
                no limitation

            TSO cost : 0
            Market cost : 0
        =#
        @testset "no_problem" begin

            context = create_instance(TS,
                                        [10., 25.], #load1
                                        [30., 35.], #load2
                                        [40., 60.], #wind_1_1
                                    35.)

            tso = PSCOPF.TSOBilevel()
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            test_ts1_values(result, TS)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[2],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            #Market EOD constraints are OK
            @test value(result.lower.limitable_model.p_global_capping[TS[2],"S1"]) < 1e-09
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            # no limitation
            @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[2],"S1"]) > 60. - 1e-09
            @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[2],"S1"]) < 1e-09
            @test value(result.upper.limitable_model.p_capping["wind_1_1",TS[2],"S1"]) < 1e-09

            #costs
            @test value(PSCOPF.get_upper_obj_expr(result)) < 1e-09
            @test value(PSCOPF.get_lower_obj_expr(result)) < 1e-09
        end

        #=
        TS: [11h, 11h15]
        S: [S1]
        TS1       TS2
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 40    S1: 60    |----------------------|
                            |         35           |
                            |                      |
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 10    S1: 10    |                      | S1: 30   S1: 30

        At 11h15 :
            Available limitable production : 60
            Demand : 40

            Market :
                The market needs to satisfy EOD constraint => only option is to cap 20MW of limitable production.
                    e = p_capping = 20.
                    lol = p_global_loss_of_load = 0.

            The TSO locates the capping to avoid RSO constraints (even if there is no risk of violating RSO here)
                p_capping[wind_1_1]=20, plim[wind1]=40, p_injected[wind_1_1]=40

            The TSO takes a limitation action:
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            plim[wind1] = 40 (due to locating capping)

            TSO cost : 0.001 (one limitation)
            Market cost : 20 (capping)

            NOTE: here, there is only one limitable,
            normally the TSO should not have to limit the production of wind_1_1
            cause the market does not have other options
            FIXME? : TSO limits for EOD reasons.



        In TS2, we have excess production
        => Market should cap the production.
        => limitation of the prod
        Note: the TSO does a limitation even if it's due to EOD constraints
        =#
        @testset "eod_overproduction_requires_market_capping_and_tso_limitation_in_ts2" begin

            context = create_instance(TS,
                                        [10., 10.], #load1
                                        [30., 30.], #load2
                                        [40., 60.], #wind_1_1
                                    35.)

            tso = PSCOPF.TSOBilevel()
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            test_ts1_values(result, TS)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            #TSO RSO constraints are OK => no minimum capping or cut conso
            @test value(result.upper.limitable_model.p_global_capping[TS[2],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            # TSO needs to limit because it is him who locates capping
            @test 40. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[2],"S1"])
            @test 1. ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[2],"S1"])
            @test 20. ≈ ( value(result.upper.limitable_model.p_capping["wind_1_1",TS[2],"S1"]))

            # EOD problem in TS2 : lots of prod
            @test 20. ≈ value(result.lower.limitable_model.p_global_capping[TS[2],"S1"])
            @test value(result.lower.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            #costs
            @test (1*tso.configs.TSO_LIMIT_PENALTY) ≈ value(PSCOPF.get_upper_obj_expr(result)) #one limitation
            @test (20. *tso.configs.MARKET_CAPPING_COST) ≈ value(PSCOPF.get_lower_obj_expr(result))

        end

        #=
        TS: [11h, 11h15]
        S: [S1]
        TS1       TS2
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 40    S1: 40    |----------------------|
                            |         35           |
                            |                      |
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 10    S1: 30    |                      | S1: 30   S1: 30

        At 11h15:

            Available limitable production : 40
            Demand : 60

            Market :
                The market needs to satisfy EOD => reduces the consumption to 40
                    e = p_capping = 0.
                    lol = p_loss_of_load = 20.

            The TSO locates the LoL to avoid RSO constraints
                Here, there is no risk of violating the RSO constraint => any option will do:
                    p_loss_of_load[bus_1]+p_loss_of_load[bus_2] = 20

            The TSO takes no further actions
            e_min = p_global_capping = 0.
            lol_min = p_global_loss_of_load = 0.
            no limitation

            TSO cost : 0
            Market cost : 20*1e5 (capping)
        =#
        @testset "eod_overload_requires_market_lol_in_ts2" begin

            context = create_instance(TS,
                                        [10., 30.], #load1
                                        [30., 30.], #load2
                                        [40., 40.], #wind_1_1
                                    35.)

            tso = PSCOPF.TSOBilevel()
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            test_ts1_values(result, TS)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

            #TSO RSO constraints are OK
            @test value(result.upper.limitable_model.p_global_capping[TS[2],"S1"]) < 1e-09
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            #no limitation
            @test value(result.upper.limitable_model.p_limit["wind_1_1",TS[2],"S1"]) > 40. - 1e-09
            @test value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[2],"S1"]) < 1e-09
            @test value(result.upper.limitable_model.p_capping["wind_1_1",TS[2],"S1"]) < 1e-09

            # EOD problem in TS2 : lot of conso
            @test value(result.lower.limitable_model.p_global_capping[TS[2],"S1"]) < 1e-09
            @test 20. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[2],"S1"])

            #TSO distributes the cut conso assuring RSO
            @test 20. ≈ ( value(result.upper.lol_model.p_loss_of_load["bus_1", TS[2],"S1"])
                        + value(result.upper.lol_model.p_loss_of_load["bus_2", TS[2],"S1"]) )

            #costs
            @test (20*tso.configs.TSO_LOL_PENALTY) ≈ value(PSCOPF.get_upper_obj_expr(result))
            @test (20*tso.configs.MARKET_LOL_PENALTY) ≈ value(PSCOPF.get_lower_obj_expr(result))

        end

        #=
        TS: [11h, 11h15]
        S: [S1]
        TS1       TS2
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 40    S1: 50    |----------------------|
                            |         35           |
                            |                      |
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 10    S1: 10    |                      | S1: 30   S1: 40

        At 11h15 :
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
        =#
        @testset "rso_problem_in_ts2" begin

            context = create_instance(TS,
                                        [10., 10.], #load1
                                        [30., 40.], #load2
                                        [40., 50.], #wind_1_1
                                    35.)

            tso = PSCOPF.TSOBilevel()
            firmness = PSCOPF.compute_firmness(tso,
                                                ech, next_ech,
                                                TS, context)
            result = PSCOPF.run(tso, ech, firmness,
                        PSCOPF.get_target_timepoints(context),
                        context)

            test_ts1_values(result, TS)

            # Market needed to cut conso
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

            #TSO RSO constraints are violated in TS2 => need to limit flow to 35
            @test 5. ≈ value(result.upper.limitable_model.p_global_capping[TS[2],"S1"])
            @test value(result.upper.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

            #limit wind_1_1 at TS2
            @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[2],"S1"])
            @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[2],"S1"])
            @test 5. ≈ value(result.upper.limitable_model.p_capping["wind_1_1",TS[2],"S1"])

            # EOD problem in TS2 due to TSO action solved by cutting conso
            @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[2],"S1"])
            @test 5. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[2],"S1"])

            #TSO distributes the cut conso assuring RSO
            @test value(result.upper.lol_model.p_loss_of_load["bus_1", TS[2],"S1"]) < 1e-09
            @test 5. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2", TS[2],"S1"])

            #costs
            @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (5*tso.configs.TSO_CAPPING_COST
                                                                + 1*tso.configs.TSO_LIMIT_PENALTY
                                                                + 5*tso.configs.TSO_LOL_PENALTY)
            @test (5*tso.configs.MARKET_CAPPING_COST + 5*tso.configs.MARKET_LOL_PENALTY) ≈ value(PSCOPF.get_lower_obj_expr(result))
        end

    end

    #=
    TS: [11h, 11h15]
    S: [S1]
    TS1       TS2
                    bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
    S1: 60    S1: 50    |----------------------|
                        |         35           |
                        |                      |
                        |                      |
        load(bus_1)     |                      |load(bus_2)
    S1: 15    S1: 10    |                      | S1: 45   S1: 40

    At 11h15 :
        Available limitable production : 60
        Demand : 60

        Market :
            The EOD is satisfied, the market would simply use the limitable production to satisfy the demand.
            => this will cause a flow[branch_1_2] = 45 > 35  ==> violates the RSO constraint

        The TSO needs to impose some additional constraints to oblige the market to respect the RSO constraint:
            option 1 : impose a capping of 10MW, limit wind_1_1's production to 50MW
            => cost : 10 + 0.001
                This will oblige the market to cut conso by 10.
                Since the TSO locates the cut conso, he can assure cutting conso on bus2 to assure the RSO constraint
            option 2 : impose reducing the consumption of bus2 by 10MW, limit wind_1_1's production to 50MW
                bus 2 will only require 35MW which assures RSO constraints
                TSO will need to limit the prod of wind_1_1
            => cost : 10e05 + 0.001

        TSO adopts the cheaper option : option 1
            e_min = p_global_capping = 10.
            lol_min = p_global_loss_of_load = 0.
            plimit[wind_1_1] = 50

        Market
        The market is now constrained by the TSO decisions:
            e = p_capping = 10. (due to TSO decisions)
            lol = p_loss_of_load = 10. (due to EOD)

        The TSO locates the capping and LoL:
            p_capping[wind_1_1]=10
            p_loss_of_load[bus_1]=0, p_loss_of_load[bus_2]=10

        TSO cost : 10.001
        Market cost : 10 + 10*1e5

    At 11h15 :
        c.f. rso_problem_in_ts2
        TSO cost : 5.001
        Market cost : 5 + 5*1e5

    Total cost :
        TSO cost : 15.002
        Market cost : 15 + 15*1e5

    =#
    @testset "rso_problem_in_both_ts" begin

        context = create_instance(TS,
                                    [15., 10.], #load1
                                    [45., 40.], #load2
                                    [60., 50.], #wind_1_1
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

        #TS1

        #TSO RSO constraints are violated in TS1 => need to limit flow to 35
        @test 10. ≈ value(result.upper.limitable_model.p_global_capping[TS[1],"S1"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #limit wind_1_1 at TS1
        @test 50. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[1],"S1"])
        @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[1],"S1"])
        @test 10. ≈ value(result.upper.limitable_model.p_capping["wind_1_1",TS[1],"S1"])

        # EOD problem in TS1 due to TSO action solved by cutting conso
        @test 10. ≈ value(result.lower.limitable_model.p_global_capping[TS[1],"S1"])
        @test 10. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"])

        #TSO distributes the cut conso assuring RSO
        @test value(result.upper.lol_model.p_loss_of_load["bus_1", TS[1],"S1"]) < 1e-09
        @test 10. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2", TS[1],"S1"])

        #TS2

        #TSO RSO constraints are violated in TS2 => need to limit flow to 35
        @test 5. ≈ value(result.upper.limitable_model.p_global_capping[TS[2],"S1"])
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[2],"S1"]) < 1e-09

        #limit wind_1_1 at TS2
        @test 45. ≈ value(result.upper.limitable_model.p_limit["wind_1_1",TS[2],"S1"])
        @test 1 ≈ value(result.upper.limitable_model.b_is_limited["wind_1_1",TS[2],"S1"])
        @test 5. ≈ value(result.upper.limitable_model.p_capping["wind_1_1",TS[2],"S1"])

        # EOD problem in TS2 due to TSO action solved by cutting conso
        @test 5. ≈ value(result.lower.limitable_model.p_global_capping[TS[2],"S1"])
        @test 5. ≈ value(result.lower.lol_model.p_global_loss_of_load[TS[2],"S1"])

        #TSO distributes the cut conso assuring RSO
        @test value(result.upper.lol_model.p_loss_of_load["bus_1", TS[2],"S1"]) < 1e-09
        @test 5. ≈ value(result.upper.lol_model.p_loss_of_load["bus_2", TS[2],"S1"])


        #costs
        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (
                    10*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY + 10*tso.configs.TSO_LOL_PENALTY #TS1
                    + 5*tso.configs.TSO_CAPPING_COST + 1*tso.configs.TSO_LIMIT_PENALTY  + 5*tso.configs.TSO_LOL_PENALTY #TS2
                    )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (
                    10*tso.configs.MARKET_CAPPING_COST + 10*tso.configs.MARKET_LOL_PENALTY #TS1
                    + 5*tso.configs.MARKET_CAPPING_COST + 5*tso.configs.MARKET_LOL_PENALTY #TS1
                    )

    end

end
