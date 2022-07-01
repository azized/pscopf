using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_balance_market" begin

    #=
    ECH = 10h
    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 50            |----------------------|
      S2: 40            |         35           |
                        |                      |
                        |                      |
    (pilotable) prod_1_1|                      |(pilotable) prod_2_1
    Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
    Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
    ON                  |                      | ON
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 100           |                      | S1: 40
      S2: 100           |                      | S2: 40
    =#
    function create_instance(ech, next_ech, TS;
                            limit_1::Union{Missing,Float64}=missing,
                            impos_1_s1::Union{Missing,Tuple{Float64,Float64}}=missing, impos_1_s2::Union{Missing,Tuple{Float64,Float64}}=missing,
                            impos_2_s1::Union{Missing,Tuple{Float64,Float64}}=missing, impos_2_s2::Union{Missing,Tuple{Float64,Float64}}=missing,
                            )
        network = PSCOPFFixtures.network_2buses()
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                10., 100.,
                                                450., 10.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                10., 100.,
                                                800., 15.,
                                                Dates.Second(0), Dates.Second(0))
        # initial generators state
        generators_init_state = SortedDict(
                        "prod_1_1" => PSCOPF.ON,
                        "prod_2_1" => PSCOPF.ON
                    )
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], "S1", 50.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", TS[1], "S2", 40.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S2", 100.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", 40.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S2", 40.)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing,
                                        "debug_balance")

        tso_actions_l = PSCOPF.get_tso_actions(context)
        if !ismissing(limit_1)
            PSCOPF.set_limitation_definitive_value!(tso_actions_l, "wind_1_1", TS[1], limit_1)
        end
        if !ismissing(impos_1_s1)
            PSCOPF.set_imposition_value!(tso_actions_l, "prod_1_1", TS[1], "S1", impos_1_s1[1], impos_1_s1[2])
        end
        if !ismissing(impos_1_s2)
            PSCOPF.set_imposition_value!(tso_actions_l, "prod_1_1", TS[1], "S2", impos_1_s2[1], impos_1_s2[2])
        end
        if !ismissing(impos_2_s1)
            PSCOPF.set_imposition_value!(tso_actions_l, "prod_2_1", TS[1], "S1", impos_2_s1[1], impos_2_s1[2])
        end
        if !ismissing(impos_2_s2)
            PSCOPF.set_imposition_value!(tso_actions_l, "prod_2_1", TS[1], "S2", impos_2_s2[1], impos_2_s2[2])
        end

        market = PSCOPF.BalanceMarket()
        firmness = PSCOPF.compute_firmness(market, ech, next_ech, TS, context)
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        return context, result
    end

    ech = DateTime("2015-01-01T10:00:00")
    TS = [DateTime("2015-01-01T11:00:00")]

    #=
    No limit, No impositions

    ECH = 10h
    TS: [11h]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=100    |                      |
    Csta=0, Cprop=1     |                      |
    S1: 50              |----------------------|
    S2: 40              |         35           |
                        |                      |
                        |                      |
    (pilotable) prod_1_1|                      |(pilotable) prod_2_1
    Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
    Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
    ON                  |                      | ON
                        |                      |
        load(bus_1)     |                      |load(bus_2)
    S1: 100             |                      | S1: 40
    S2: 100             |                      | S2: 40

    wind_1_1 will produce at maximum allowed ;
    S1:50 , S2:40
    prod_1_1 is cheaper => will be used :
    S1:90 , S2:100
    prod_2_1 : shutdown
    =#
    @testset "no_limitations_no_impositions" begin
        context, result = create_instance(ech, ech+Minute(1), TS)
        m_schedule_l = PSCOPF.get_market_schedule(context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # Limitable respects limitations
        @test 50. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
        @test 40. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
        @test value(result.limitable_model.p_global_capping[TS[1], "S1"]) < 1e-09
        @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09

        @test 90. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
        @test 100. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")
        @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1") < 1e-09
        @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2") < 1e-09
    end

    @testset "balance_market_respects_limitations" begin

        #=
        Limit wind_1_1 to 45MW => only limits S1

        ECH = 10h
        TS: [11h]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 50              |----------------------|
        S2: 40              |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
        ON                  |                      | ON
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 100             |                      | S1: 40
        S2: 100             |                      | S2: 40

        wind_1_1 will produce at maximum allowed ;
        S1:45 , S2:40
        prod_1_1 is cheaper => will be used :
        S1:95 , S2:100
        prod_2_1 : shutdown
        =#
        @testset "limit_affects_one_scenario" begin
            context, result = create_instance(ech, ech+Minute(1), TS, limit_1=45.)
            m_schedule_l = PSCOPF.get_market_schedule(context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            # Limitable respects limitations
            @test 45. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
            @test 40. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
            @test 5. ≈ value(result.limitable_model.p_global_capping[TS[1], "S1"])
            @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09

            @test 95. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
            @test 100. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")
            @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1") < 1e-09
            @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2") < 1e-09
        end

        #=
        Limit wind_1_1 to 35MW => limits S1 and S2

        ECH = 10h
        TS: [11h]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 50              |----------------------|
        S2: 40              |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
        ON                  |                      | ON
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 100             |                      | S1: 40
        S2: 100             |                      | S2: 40

        wind_1_1 will produce at maximum allowed ;
        S1:35 , S2:35
        prod_1_1 is cheaper => want to use it at 100MW but prod_2_1 has a minimum capacity of 10MW :
        S1:95 , S2:95
        prod_2_1 :
        S1:10 , S2:10
        =#
        @testset "limit_affects_both_scenarios" begin
            context, result = create_instance(ech, ech+Minute(1), TS, limit_1=35.)
            m_schedule_l = PSCOPF.get_market_schedule(context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            # Limitable respects limitations
            @test 35. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
            @test 35. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
            @test 15. ≈ value(result.limitable_model.p_global_capping[TS[1], "S1"]) #35MW /50
            @test 5. ≈ value(result.limitable_model.p_global_capping[TS[1], "S2"]) #35MW /40

            @test 95. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
            @test 95. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")
            @test 10. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1")
            @test 10. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2")
        end

    end

    @testset "balance_market_respects_impositions" begin

        #=
        impose on prod_1_1 to be in:
        S1 : [10, 85]
        S2 : [10, 95]
        The market can no longer choose to shutdown prod_1_1

        ECH = 10h
        TS: [11h]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 50              |----------------------|
        S2: 40              |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
        ON                  |                      | ON
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 100             |                      | S1: 40
        S2: 100             |                      | S2: 40

        wind_1_1 will produce at maximum allowed ;
            S1:50 , S2:40
        prod_1_1 is cheaper => will be used as allowed by impositions:
            S1:80 (not 85 cause need to use prod_2_1 and prod_2_1 has a min capacity of 10),
            S2:90 (not 95 cause need to use prod_2_1 and prod_2_1 has a min capacity of 10),
        prod_2_1 at p_min
            S1:10 , S2:10
        =#
        @testset "imposition_on_prod_1_1" begin
            impos_1_s1_l = (10., 85.)
            impos_1_s2_l = (10., 95.)

            context, result = create_instance(ech, ech+Minute(1), TS,
                                            impos_1_s1=impos_1_s1_l,
                                            impos_1_s2=impos_1_s2_l,
                                            )
            m_schedule_l = PSCOPF.get_market_schedule(context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            # Pilotable respects impositions
            @test 80. ≈ value(result.pilotable_model.p_injected["prod_1_1",TS[1], "S1"])
            @test 80. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1"), impos_1_s1_l[1], impos_1_s1_l[2])
            @test 90. ≈ value(result.pilotable_model.p_injected["prod_1_1",TS[1], "S2"])
            @test 90. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2"), impos_1_s2_l[1], impos_1_s2_l[2])

            @test 10. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1")
            @test 10. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2")

            # Limitable is used at maximum allowed
            @test 50. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
            @test 40. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
            @test value(result.limitable_model.p_global_capping[TS[1], "S1"]) < 1e-09
            @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09
        end

        #=
        impose on prod_1_1 to be in:
        S1 : [0, 75]
        S2 : [0, 80]
        The market can choose to shutdown prod_1_1 or not
        Here it will use it cause it's cheaper than prod_2_1

        ECH = 10h
        TS: [11h]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 50              |----------------------|
        S2: 40              |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
        ON                  |                      | ON
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 100             |                      | S1: 40
        S2: 100             |                      | S2: 40

        wind_1_1 will produce at maximum allowed ;
            S1:50 , S2:40
        prod_1_1 is cheaper => will be used as allowed by impositions:
            S1:75
            S2:80
        prod_2_1 at p_min
            S1:15 , S2:20
        =#
        @testset "different_impositions_by_scenario" begin
            impos_1_s1_l=(0., 75.)
            impos_1_s2_l=(0., 80.)

            context, result = create_instance(ech, ech+Minute(1), TS,
                                            impos_1_s1=impos_1_s1_l,
                                            impos_1_s2=impos_1_s2_l,
                                            )
            m_schedule_l = PSCOPF.get_market_schedule(context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            # Pilotable respects impositions
            @test 75. ≈ value(result.pilotable_model.p_injected["prod_1_1",TS[1], "S1"])
            @test 75. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1"), impos_1_s1_l[1], impos_1_s1_l[2])
            @test 80. ≈ value(result.pilotable_model.p_injected["prod_1_1",TS[1], "S2"])
            @test 80. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2"), impos_1_s2_l[1], impos_1_s2_l[2])

            @test 15. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1")
            @test 20. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2")

            # Limitable is used at maximum allowed
            @test 50. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
            @test 40. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
            @test value(result.limitable_model.p_global_capping[TS[1], "S1"]) < 1e-09
            @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09
        end

        #=
        impose on prod_2_1 to be in:
        S1 : [0, 75]
        S2 : [0, 80]
        The market can choose to shutdown prod_2_1 or not
        Here it will shut it down cause it's cheaper to use prod_1_1

        ECH = 10h
        TS: [11h]
        S: [S1,S2]
                            bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=100    |                      |
        Csta=0, Cprop=1     |                      |
        S1: 50              |----------------------|
        S2: 40              |         35           |
                            |                      |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=10, Pmax=100   |                      | Pmin=10, Pmax=100
        Csta=450, Cprop=10  |                      | Csta=800, Cprop=15
        ON                  |                      | ON
                            |                      |
            load(bus_1)     |                      |load(bus_2)
        S1: 100             |                      | S1: 40
        S2: 100             |                      | S2: 40

        wind_1_1 will produce at maximum allowed ;
            S1:50 , S2:40
        prod_1_1 is cheaper => will be used as allowed by impositions:
            S1:90
            S2:100
        prod_2_1 at p_min
            S1:0 , S2:0
        =#
        @testset "imposition_levels_may_leave_commitment_decision_to_market" begin
            impos_2_s1_l=(0., 75.)
            impos_2_s2_l=(0., 80.)

            context, result = create_instance(ech, ech+Minute(1), TS,
                                            impos_2_s1=impos_2_s1_l,
                                            impos_2_s2=impos_2_s2_l,
                                            )
            m_schedule_l = PSCOPF.get_market_schedule(context)

            # Solution is optimal
            @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

            # Pilotable respects impositions
            @test 90. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S1")
            @test 100. ≈ PSCOPF.get_prod_value(m_schedule_l, "prod_1_1", TS[1], "S2")

            @test value(result.pilotable_model.p_injected["prod_2_1",TS[1], "S1"]) < 1e-09
            @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1") < 1e-09
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S1"), impos_2_s1_l[1], impos_2_s1_l[2])
            @test value(result.pilotable_model.p_injected["prod_2_1",TS[1], "S2"])  < 1e-09
            @test PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2")  < 1e-09
            @test PSCOPF.in_bounds(PSCOPF.get_prod_value(m_schedule_l, "prod_2_1", TS[1], "S2"), impos_2_s2_l[1], impos_2_s2_l[2])

            # Limitable is used at maximum allowed
            @test 50. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S1")
            @test 40. ≈ PSCOPF.get_prod_value(m_schedule_l, "wind_1_1", TS[1], "S2")
            @test value(result.limitable_model.p_global_capping[TS[1], "S1"]) < 1e-09
            @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09
        end

    end

end
