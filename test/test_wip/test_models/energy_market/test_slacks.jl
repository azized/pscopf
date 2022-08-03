using PSCOPF

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_energy_market_slacks" begin

    #=
    If available power for a limitable exceeds consumption,
      excess power is capped.
    The capping variable is global (delocalised in space but has
      a per timestep and per scenario value)
    The power capped is distributed in the schedule (not visible here cause only one limitable unit).

    TS: [11h]
    S: [S1]
                      bus 1
                        |
    (limitable) wind_1_1|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |     load_1
      S1: 20            | S1: 15
      S2: 25            | S2: 25
                        |
    =#
    @testset "energy_market_capping_limitable_power" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 20.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 15.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict{String, SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness} }(),
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE))
                    )
        # initial generators state : No need because all pmin=0 => ON by default
        generators_init_state = SortedDict{String, PSCOPF.GeneratorState}()

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        # Limitable produces to the available level
        @test 15. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S2")

        # Capping is not localised in market model
        # Limitable was capped when prod > load (ie. S1):
        @test 5. ≈ value(result.limitable_model.p_global_capping[TS[1], "S1"])
        # Limitable was not capped when prod <= load (ie. S2):
        @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09

        # Capping is localised in the market schedule
        @test 5. ≈ PSCOPF.get_capping(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test PSCOPF.get_capping(context.market_schedule, "wind_1_1", TS[1], "S2") < 1e-09

        # No penalty but we do pay for capped power
        #(otherwise we need a per unit variable for capping to localise and reduce the unpaid costs)
        @test value(result.objective_model.penalty) < 1e-09
        @test value(result.objective_model.start_cost) < 1e-09
        @test value(result.objective_model.prop_cost) < 1e-09 # no cost for limitables
    end


    #=
    It is possible to cut consumption.
    This variable should can be used in two cases :
    1- due to EOD constraints : we don't have enough production capacity (illustrated in S1)
    2- due to Pmin constraints : we don't have enough demand to start a unit (illustrated in S2)

    TS: [11h]
    S: [S1]
                        bus 1
                         |
    (pilotable) prod_1_1 |    load_1
    Pmin=20, Pmax=100    |  S1: 150
    Csta=1000, Cprop=1   |  S2: 15
                         |  S3: 25
    =#
    @testset "energy_market_cuts_consumption_due_to_pmin" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 100.,
                                                1000., 1.,
                                                Dates.Second(0), Dates.Second(0))
        prod_1_1 = PSCOPF.safeget_generator(network, "prod_1_1")
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 150.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 15.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S3", 25.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),),
                    SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE))
                    )
        # initial generators state :
        generators_init_state = SortedDict("prod_1_1" => PSCOPF.OFF)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing, "DBGUG")
        market = PSCOPF.EnergyMarket()
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # indicates using slacks for feasibility
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK
        # S1 : prod_capacity < load => cannot satisfy demand
        @test 100. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1")
        @test 50. ≈ value(result.lol_model.p_global_loss_of_load[TS[1], "S1"])
        # S2 : load < pmin => cannot start the unit for such load
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2") < 1e-09
        @test value(result.limitable_model.p_global_capping[TS[1], "S2"]) < 1e-09
        @test 15. ≈ value(result.lol_model.p_global_loss_of_load[TS[1], "S2"])
        # S3 : works fine
        @test 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S3")
        @test value(result.lol_model.p_global_loss_of_load[TS[1], "S3"]) < 1e-09
        # penalize cutting consumption
        lol_cost = market.configs.loss_of_load_penalty
        @test (50. * lol_cost + 15. * lol_cost + 0. ) ≈ value(result.objective_model.penalty)
        @test (2 * PSCOPF.get_start_cost(prod_1_1)) ≈ value(result.objective_model.start_cost)
        @test ((100. + 0. + 25.)* PSCOPF.get_prop_cost(prod_1_1) ) ≈ value(result.objective_model.prop_cost)

        #cut conso is localized in the ranscripted schedule:
        @test 50. ≈ PSCOPF.get_loss_of_load(context.market_schedule, "bus_1", TS[1], "S1")
        @test 15. ≈ PSCOPF.get_loss_of_load(context.market_schedule, "bus_1", TS[1], "S2")
        @test PSCOPF.get_loss_of_load(context.market_schedule, "bus_1", TS[1], "S3") < 1e-09
    end

    #=
    Cutting consumption is done through penalization.
    Careful for the chosen penalty cost with respect to production costs
    Here, the unit starting cost is higher than the penalty

    TS: [11h]
    S: [S1]
                        bus 1
                         |
    (pilotable) prod_1_1 |    load_1
    Pmin=20, Pmax=100    |  S1: 25
    Csta=100k, Cprop=1   |

    Desired solution :
        prod_1_1 at 25MW
    Retrieved solution with low penalization (1e3):
        prod_1_1 at 0
        loss_of_load : 25MW
    =#
    @testset "energy_market_careful_for_cuts_consumption_penalty" begin
        penalty = 1e3

        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 100.,
                                                100000., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 25.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),),
                    SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE))
                    )
        # initial generators state :
        generators_init_state = SortedDict("prod_1_1" => PSCOPF.OFF)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        market.configs.loss_of_load_penalty = penalty
        @test_throws AssertionError ( result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context) )
        # PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Desired Solution
        # @test !( PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL )
        # @test !( 25. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") )
        # @test !( value(result.lol_model.p_global_loss_of_load[TS[1], "S1"]) < 1e-09 )
        # @test !( value(result.objective_model.penalty) < 1e-09 )

        # retrieved solution
        # @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK
        # @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        # @test 25. ≈ value(result.lol_model.p_global_loss_of_load[TS[1], "S1"])
        # @test (25. * 1e3) ≈ value(result.objective_model.penalty)

    end


    #=
    In S1,
    demand = 15, wind_1_1 provides 10 => residual demand of 5MW but pmin(prod_1_1)=20
    => start prod_1_1 at 20MW exceeds the demand
    => use the 10MW from wind_1_1 and cut the remaining 5MW (ie LoL)
    In S2,
    demand = 25, wind_1_1 provides 10 => residual demand of 15MW but pmin(prod_1_1)=20
    => start prod_1_1 at 20MW
    => use 5MW from wind_1_1 and cap 5MW

    TS: [11h]
    S: [S1]
                        bus 1
                         |
    (limitable) wind_1_1 |    load_1
    Pmin=0, Pmax=100     |  S1: 15
    Csta=0, Cprop=0.     |  S2: 25
         S1 : 10         |
         S1 : 10         |
                         |
    (pilotable) prod_1_1 |
    Pmin=20, Pmax=100    |
    Csta=100k, Cprop=100 |
    =#
    @testset "energy_market_capping_limitables_due_to_pilotable_pmin" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 100.,
                                                1000., 20.,
                                                Dates.Second(0), Dates.Second(0))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 10.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S2", 10.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 15.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict("prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),),
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                               "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE))
                    )
        # initial generators state :
        generators_init_state = SortedDict("prod_1_1" => PSCOPF.OFF)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # indicates using slacks for feasibility because of scenario S1
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        # In S1 : Load=15, wind provides 10 => still missing 5 but pmin=20
        # => reduce consumption by 5
        @test 10. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S1") < 1e-09
        @test value(result.limitable_model.p_global_capping[TS[1], "S1"]) < 1e-09
        @test 5. ≈ value(result.lol_model.p_global_loss_of_load[TS[1], "S1"])

        # In S2 : Load=25, wind provides 10 => still missing 15 but pmin=20
        # pilotable produces 20 => 5 extra prod (20+10 - 25) => reduce wind by 5.
        @test 5. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S2")
        @test 20. ≈ PSCOPF.get_prod_value(context.market_schedule, "prod_1_1", TS[1], "S2")
        @test 5. ≈ value(result.limitable_model.p_global_capping[TS[1], "S2"])
        @test value(result.lol_model.p_global_loss_of_load[TS[1], "S2"]) < 1e-09
    end

    #=
    For now,
     capped power is distributed evenly (in terms of percentage of available power)
     among the limitable generators
    TODO? base ratio on prop_cost ? in that case, careful for null and negative prop_cost

    TS: [11h]
    S: [S1]
                      bus 1
                        |
    (limitable) wind_1_1|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |     load_1
      S1: 50            | S1: 15
                        |
    (limitable) wind_1_2|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=2     |
      S1: 10            |
                        |
    (limitable) wind_1_3|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |
      S1: 0             |
    =#
    @testset "energy_market_capping_distribution_in_schedule" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_2", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 2.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_3", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 50.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_2", DateTime("2015-01-01T11:00:00"), "S1", 10.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_3", DateTime("2015-01-01T11:00:00"), "S1", 0.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 15.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict{String, SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness} }(),
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                                "wind_1_2" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE),
                                "wind_1_3" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE)),
                    )
        # initial generators state : No need because all pmin=0 => ON by default
        generators_init_state = SortedDict{String, PSCOPF.GeneratorState}()

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL
        # Limitable produces to the available level
        @test (0.25 * 50.) ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")
        @test (0.25 * 10.) ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_2", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.market_schedule, "wind_1_3", TS[1], "S1") <= 1e-09

        # Total capped power
        @test 45. ≈ value(result.limitable_model.p_global_capping[TS[1], "S1"])

        # Capping distributed
        @test 0.75 ≈ 45 / 60 #75% of limitable power will be capped
        @test (0.75 * 50)  ≈ context.market_schedule.capping["wind_1_1", TS[1], "S1"]
        @test (0.75 * 10.) ≈ context.market_schedule.capping["wind_1_2", TS[1], "S1"]
        @test (0.75 * 0.)  ≈ context.market_schedule.capping["wind_1_3", TS[1], "S1"]
    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=200    |                      |
    Csta=0, Cprop=1     |                      |
      S1: 160           |----------------------|
                        |         5            |
                        |                      |
    load                |                      |load
      S1: 180           |                      |  S1: 20

    available prod : 160
    total load : 200
    => need to cut 40 => 20% of load
    cut 20% on each bus :
        bus1: 20% of 180 => 36
        bus2: 20% of 20  => 4
    =#
    @testset "energy_market_lol_distribution_by_bus_in_schedule" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        # Buses
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        PSCOPF.Networks.add_new_bus!(network, "bus_2")
        # Branches
        PSCOPF.Networks.add_new_branch!(network, "branch_1_2", 5.);
        # PTDF
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_1", 0.5)
        PSCOPF.Networks.add_ptdf_elt!(network, "branch_1_2", "bus_2", -0.5)
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 200.,
                                                0., 1.,
                                                Dates.Second(0), Dates.Second(0))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 160.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 180.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", DateTime("2015-01-01T11:00:00"), "S1", 20.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict{String, SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness} }(),
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE)),
                    )
        # initial generators state : No need because all pmin=0 => ON by default
        generators_init_state = SortedDict{String, PSCOPF.GeneratorState}()

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        market = PSCOPF.EnergyMarket()
        result = PSCOPF.run(market, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_market_schedule!(context, ech, result, firmness, market)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK
        # Limitable produces to the available level
        @test 160. ≈ PSCOPF.get_prod_value(context.market_schedule, "wind_1_1", TS[1], "S1")

        # Total cut conso
        @test 40. ≈ value(result.lol_model.p_global_loss_of_load[TS[1], "S1"])

        # cut conso distributed on buses
        @test 0.2 ≈ 40 / 200 # 40:loss_of_load, 200:total load
        @test (0.2 * 180) ≈ context.market_schedule.loss_of_load_by_bus["bus_1", TS[1], "S1"]
        @test (0.2 * 20.) ≈ context.market_schedule.loss_of_load_by_bus["bus_2", TS[1], "S1"]
    end

    #TODO : illustrate ptdf effect
    #TODO : test capping and load_shedding transcription in schedule for marketAtFO

end
