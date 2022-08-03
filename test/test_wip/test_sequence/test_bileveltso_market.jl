using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bileveltso_market" begin

    function create_instance(ECH, S, TS,
                            fo::Minute,
                            load_1::Vector{Vector{Float64}}, load_2::Vector{Vector{Float64}}, wind_1::Vector{Vector{Float64}},
                            preceding_decision_prod_1::Vector{Vector{PSCOPF.GeneratorState}},
                            preceding_decision_prod_2::Vector{Vector{PSCOPF.GeneratorState}},
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 200.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 200.,
                                                100., 10.,
                                                Dates.Second(4*60*60), Dates.Second(15*60))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                20., 200.,
                                                500., 50.,
                                                Dates.Second(1*60*60), Dates.Second(15*60))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()

        for (s_index, s) in enumerate(S)
            for (ts_index, ts) in enumerate(TS)
                PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_1",  ts, s, load_1[s_index][ts_index]  )
                PSCOPF.add_uncertainty!(uncertainties, ECH[1], "bus_2",  ts, s, load_2[s_index][ts_index])
                PSCOPF.add_uncertainty!(uncertainties, ECH[1], "wind_1_1", ts, s, wind_1[s_index][ts_index])
            end
        end

        gen_initial_state = SortedDict{String,PSCOPF.GeneratorState}(
                                "prod_1_1" => PSCOPF.OFF,
                                "prod_2_1" => PSCOPF.OFF,
                            )

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.ManagementMode("test_mode", fo),
                                    gen_initial_state,
                                    uncertainties, nothing, logs)

        # c.f. tso_bilevel_does_not_start_units_for_eod_reasons_from_missing
        # important in context.market_schedule for TSOBilevel, because tso bilevel won't start units for EOD reasons
        for (s_index, s) in enumerate(S)
            for (ts_index, ts) in enumerate(TS)
                PSCOPF.set_commitment_value!(context.market_schedule, "prod_1_1", ts, s, preceding_decision_prod_1[s_index][ts_index])
                PSCOPF.set_commitment_value!(context.tso_schedule, "prod_1_1", ts, s, preceding_decision_prod_1[s_index][ts_index])
                PSCOPF.set_commitment_value!(context.market_schedule, "prod_2_1", ts, s, preceding_decision_prod_2[s_index][ts_index])
                PSCOPF.set_commitment_value!(context.tso_schedule, "prod_2_1", ts, s, preceding_decision_prod_2[s_index][ts_index])
            end
        end

        PSCOPF.set_horizon_timepoints(context, ECH)
        PSCOPF.check(context)

        return context
    end




################################################################################################
#           TEST CASE
################################################################################################

    #=
        TS: [11h]
        S: [S1]
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=20, Pmax=200   |                      |
        Csta=100, Cprop=10  |                      |
        S1: 85              |----------------------|
                            |         35           |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=20, Pmax=200   |                      | Pmin=20, Pmax=200
        Csta=500, Cprop=50  |                      | Csta=50k, Cprop=50
        OFF->ON             |                      | OFF->OFF
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 50              |                      | S1: 70
    =#
    @testset "test_two_modes_on_a_single_ech_a_single_ts" begin
        ECH = [DateTime("2015-01-01T07:00:00"), DateTime("2015-01-01T10:00:00")]
        S = ["S1"]
        TS = [DateTime("2015-01-01T11:00:00")]
        future_ech = ECH[2]

        #[[S1_TS1, S1_TS_2, ...], [S2_TS1, S2_TS_2, ...]]
        load_1 = [[50.]]
        load_2 = [[70.]]
        wind_1 = [[85.]]
        test2_preceding_decision_prod_1 = [[PSCOPF.ON]]
        test2_preceding_decision_prod_2 = [[PSCOPF.OFF]]

        #=
        TS: [11h]
        S: [S1]
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=20, Pmax=200   |                      |
        Csta=100, Cprop=10  |                      |
        S1: 85              |----------------------|
                            |         35           |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=20, Pmax=200   |                      | Pmin=20, Pmax=200
        Csta=500, Cprop=50  |                      | Csta=50k, Cprop=50
        OFF->ON             |                      | OFF->OFF
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 50              |                      | S1: 70

        MARKET:
        demand = 120
        market chooses economically,
        It uses the limitable prod => P_wind_1_1 = 85.
        remaining demand : 120-85 = 35MW
        Then the cheaper unit (start and prop cost) => prod_1_1 => P_prod_1_1 = 35.

        This causes a flow of 70 (>35.) on branch_1_2, the TSO needs to react

        TSO:
        The tso solves the RSO constraint while trying to stick to the market levels

        option 1 :
            -35 on prod_1_1, +35 on prod_2_1
            => Delta = 70 => cost_1 = 70.
            Pwind_1_1 = 85, P_prod_2_1 = 35
            cost = 85 + 35. * 50 + 50000 = 51835
        option 2 :
            -20 on wind_1_1, -15 on prod_1_1, +35 on prod_2_1
            => Delta = 70 => cost_1 = 70. + 1e-3 (cause one limitation)
            cost_step1 = 70.001 => option 1 is better
            Pwind_1_1 = 65, P_prod_2_1 = 20, P_prod_2_1 = 35
            cost = 65 + 20. * 10 + 35. * 50 + 50000 = 52015

        Market Decision :
            P_wind_1_1 = 85
            P_prod_1_1 = 35
            P_prod_2_1 = 0
        TSO reaction :
            P_wind_1_1 = 85 (unchanged)
            P_prod_1_1 = 0 (-35)
            P_prod_2_1 = 35 (+35)

        =#
        @testset "energymarket_tso" begin

            context = create_instance(ECH, S, TS, Dates.Minute(1),
                                    load_1,
                                    load_2,
                                    wind_1,
                                    test2_preceding_decision_prod_1,
                                    test2_preceding_decision_prod_2,
                                    35.)

            @testset "energymarket_before_tso" begin
                PSCOPF.run_step!(context, PSCOPF.EnergyMarket(), ECH[1], ECH[2])

                @test 85. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
                @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
                @test PSCOPF.get_prod_value(context.market_schedule, "prod_2_1", TS[1], "S1") < 1e-09
            end

            @testset "tso" begin
                PSCOPF.run_step!(context, PSCOPF.TSOOutFO(), ECH[1], ECH[2])

                # market schedule is not changed
                @test 85. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
                @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
                @test PSCOPF.get_prod_value(context.market_schedule, "prod_2_1", TS[1], "S1") < 1e-09

                # TSO reaction
                @test 85. ≈ PSCOPF.get_prod_value(context.tso_schedule, "wind_1_1", TS[1], "S1")
                @test PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S1") < 1e-09
                @test 35. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_2_1", TS[1], "S1")

                #no limitation #Note : even if there was, they would only appear for units that reached the DP
                @test isempty(PSCOPF.get_limitations(context.tso_actions))
                #impositions only appear at DP
                @test isempty(PSCOPF.get_impositions(context.tso_actions))
                @test ismissing(PSCOPF.get_imposition(context.tso_actions, "prod_1_1", TS[1], "S1"))
                @test ismissing(PSCOPF.get_imposition(context.tso_actions, "prod_2_1", TS[1], "S1"))

            end

        end

        #=
        TS: [11h]
        S: [S1]
                        bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=20, Pmax=200   |                      |
        Csta=100, Cprop=10  |                      |
        S1: 85              |----------------------|
                            |         35           |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=20, Pmax=200   |                      | Pmin=20, Pmax=200
        Csta=500, Cprop=50  |                      | Csta=50k, Cprop=50
        OFF->ON             |                      | OFF->OFF
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        S1: 50              |                      | S1: 70

        MARKET:
        demand = 120
        market chooses economically,
        It uses the limitable prod => P_wind_1_1 = 85.
        remaining demand : 120-85 = 35MW
        Then the cheaper unit => prod_1_1 => P_prod_1_1 = 35.
        This would cause a flow of 70 (>35.) on branch_1_2, the TSO needs to react

        TSO:
        The TSO needs to anticipate to prevent the market from violating the RSO constraint
        => needs to limit the flow on branch_1_2 to 35.

        option 1 : TSO cost=15.
            produce at least 35.MW on bus 2
            => P_prod_2_1 \in 35-200 => (unit was off) cost=35
            to avoid any impositions costs, P_prod_1_1 will stay at [20, 200] (imposition cost would be 180)
            Since P_prod_1_1 will produce 20MW but bus 1's production should not exceed 85MW
            TSO will also limit wind_1_1 to 65MW

            The market will need to cap 20MW on wind_1_1

        option 2 : cost=180.
            produce at most 85MW on bus 1
            => P_prod_1_1 is OFF => cost= (200-0) + (0-20) = 180
            to avoid any impositions costs, P_prod_1_1 will stay at [20, 200]

        TSO decisions : cost = 15
            P_wind_1_1 : <= 65 (due to EOD not to TSO restrictions)
            P_prod_1_1 : 20 - 200
            P_prod_2_1 : 35 - 200
        Market decisions : cost = 65 + 20*10 + 35*50 = 2015 (ignoring start cost)
            P_wind_1_1 = 65
            P_prod_1_1 = 20
            P_prod_2_1 = 35

        FIXME? we prefer ?
        TSO decisions : cost = currently 180 but should be 0 ?
            P_wind_1_1 : no limit
            P_prod_1_1 : 0 - 0 i.e. OFF
            P_prod_2_1 : 20 - 200
        Market decisions : cost = 85 + 0 + 35*50 = 1835 (ignoring start cost)
            P_wind_1_1 = 85
            P_prod_1_1 = 0.
            P_prod_2_1 = 35
        =#
        @testset "tsobilevel_balancemarket" begin

            context = create_instance(ECH, S, TS, Dates.Minute(1),
                                    load_1,
                                    load_2,
                                    wind_1,
                                    test2_preceding_decision_prod_1,
                                    test2_preceding_decision_prod_2,
                                    35.)

            @testset "tso_bilevel" begin
                PSCOPF.run_step!(context, PSCOPF.TSOBilevel(), ECH[1], ECH[2])

                @test 65. ≈ PSCOPF.get_prod_value(context.tso_schedule, "wind_1_1", TS[1], "S1")
                @test 20. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_1_1", TS[1], "S1")
                @test 35. ≈ PSCOPF.get_prod_value(context.tso_schedule, "prod_2_1", TS[1], "S1")

                #no limitation #Note : even if there was, they would only appear for units that reached the DP
                @test 65. ≈ PSCOPF.get_limitation(context.tso_actions, "wind_1_1", TS[1],"S1")
                #impositions only appear at DP
                @test !isempty(PSCOPF.get_impositions(context.tso_actions))
                @test 20. ≈ PSCOPF.get_imposition(context.tso_actions, "prod_1_1", TS[1], "S1")[1]
                @test 200. ≈ PSCOPF.get_imposition(context.tso_actions, "prod_1_1", TS[1], "S1")[2]
                @test 35. ≈ PSCOPF.get_imposition(context.tso_actions, "prod_2_1", TS[1], "S1")[1]
                @test 200. ≈ PSCOPF.get_imposition(context.tso_actions, "prod_2_1", TS[1], "S1")[2]

                @test 20. ≈ PSCOPF.get_capping(context.tso_schedule, "wind_1_1", TS[1], "S1")
                @test PSCOPF.get_loss_of_load(context.tso_schedule, "bus_1", TS[1], "S1") < 1e-09
                @test PSCOPF.get_loss_of_load(context.tso_schedule, "bus_2", TS[1], "S1") < 1e-09
            end

            @testset "balancemarket" begin
                PSCOPF.run_step!(context, PSCOPF.BalanceMarket(), ECH[1], ECH[2])

                @test 65. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
                @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
                @test 35. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_2_1", TS[1], "S1")

                @test 20. ≈ PSCOPF.get_capping(context.market_schedule, "wind_1_1", TS[1], "S1")
                @test PSCOPF.get_loss_of_load(context.market_schedule, "bus_1", TS[1], "S1") < 1e-09
                @test PSCOPF.get_loss_of_load(context.market_schedule, "bus_2", TS[1], "S1") < 1e-09
            end

        end

    end

end
