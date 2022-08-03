using PSCOPF

using Test
using JuMP
using Dates
using DataStructures

@testset verbose=true "test_tso_slacks" begin

    #=
    Currently
    Pinj = min(Plim, uncertainties)
    so Pinj cannot be fixed to a level lower than uncertainties
    ie. if limit is 25, uncertainty is 20 and demand is 15
    the expression above would cause a production of min(25,20)=20MW
    So an extra 5MW (which is not allowed in the model)
                    S1   S2
    Load            15    25
    prod            20    25
    we could have (Not what we want, injection should = min(uncertainty, Plimit)) :
    Plim         =     25
    Pinj         = [15 , 25]
    B_islim      = [1  , 0]
    Pislim_x_lim = [0  , 0]
    loss_of_load    = [0  , 0]
    We have what we want :
    Plim         =     15
    Pinj         = [15 , 15]
    B_islim      = [0  , 1]
    Pislim_x_lim = [0  , 15]
    loss_of_load    = [0  , 10]

    TS: [11h]
    S: [S1]
                      bus 1
                        |
    (limitable) wind_1_1|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |     load_1
    DP=>9h30            |
      S1: 20            | S1: 15
      S2: 25            | S2: 25
                        |
    =#
    @testset "tso_cant_cap_limitable_power_by_choosing_prod_level_after_dp" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T09:30:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(90*60), Dates.Second(90*60))
        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S1", 20.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "wind_1_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S1", 15.)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", DateTime("2015-01-01T11:00:00"), "S2", 25.)
        # firmness
        firmness = PSCOPF.Firmness(
                    SortedDict{String, SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness} }(),
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.TO_DECIDE))
                    )
        # initial generators state : No need because all pmin=0 => ON by default
        generators_init_state = SortedDict{String, PSCOPF.GeneratorState}()

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)

        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), Dates.DateTime("2015-01-01T07:00:00"), SortedDict(
            "wind_1_1" => PSCOPF.GeneratorSchedule("wind_1_1",
                SortedDict(),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>20.,"S2"=>25.)))
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK
        @test value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S1"]) < 1e-09
        @test 10. ≈ value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S2"])
        # Limitable produces to the available level
        @test 15. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S1")
        @test 15. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S2")
        # Limitable was capped when prod > load (ie. S1):
        @test -5. ≈ ( PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S1")
                    - PSCOPF.get_uncertainties(uncertainties[ech], "wind_1_1", TS[1], "S1") )
        @test -10. ≈ ( PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S2")
                    - PSCOPF.get_uncertainties(uncertainties[ech], "wind_1_1", TS[1], "S2") )
        # We do not pay for capped power
        @test 1. ≈ value(result.limitable_model.b_is_limited["wind_1_1", TS[1], "S1"])
        @test 1. ≈ value(result.limitable_model.b_is_limited["wind_1_1", TS[1], "S2"])
        @test 15. ≈ value(result.limitable_model.p_limit["wind_1_1", TS[1], "S1"])
        @test 15. ≈ value(result.limitable_model.p_limit["wind_1_1", TS[1], "S2"])
        #If it was due to TSO constraints, we would have paid for the generator used instead of limitables
        @test (10*tso.configs.loss_of_load_penalty + 2*1e-3) ≈ value(result.objective_model.penalty) # 10 loss_of_load + 2 limitations
        @test value(result.objective_model.start_cost) < 1e-09
        @test ((20. - 15.) *1. + (25. - 15.) *1. ) ≈ value(result.objective_model.prop_cost)
    end

    #=
    It is possible to cut consumption.
    This variable can be used in two cases :
    1- due to EOD constraints : we don't have enough production capacity (illustrated in S1)
    2- due to Pmin constraints : we don't have enough demand to start a unit (illustrated in S2)
    3- due to RSO constraints (c.f. tso_cuts_consumption_due_to_rso)

    TS: [11h]
    S: [S1]
                        bus 1
                         |
    (pilotable) prod_1_1 |    load_1
    Pmin=20, Pmax=100    |  S1: 150
    Csta=100k, Cprop=1   |  S2: 15
                         |  S3: 25
    =#
    @testset "tso_cuts_consumption_due_to_pmin" begin
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
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), Dates.DateTime("2015-01-01T07:00:00"), SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S1"=>PSCOPF.ON, "S2"=>PSCOPF.OFF, "S3"=>PSCOPF.ON))),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>100.,"S2"=>0.,"S3"=>25.)))
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # indicate use of slacks for feasibility
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK
        # S1 : prod_capacity < load => cannot satisfy demand
        @test 100. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "prod_1_1", TS[1], "S1")
        @test 50. ≈ value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S1"])
        # S2 : load < pmin => cannot start the unit for such load
        @test PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "prod_1_1", TS[1], "S2") < 1e-09
        @test 15. ≈ value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S2"])
        # S3 : works fine
        @test 25. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "prod_1_1", TS[1], "S3")
        @test value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S3"]) < 1e-09
        # penalize cutting consumption
        @test ((50. + 15. + 0. ) * tso.configs.loss_of_load_penalty) ≈ value(result.objective_model.penalty)
        @test (0. + 0. + 0.) ≈ value(result.objective_model.start_cost) #The market paid for starting
        @test (100. + 0. + 25. ) ≈ value(result.objective_model.prop_cost) #FIXME : pay just for the difference ? i.e. delta*prop_cost instead of injected*prop_cost
    end


    #=
    TS: [11h]
    S: [S1]
                        bus 1
                         |
    (limitable) wind_1_1 |    load_1
    Pmin=0, Pmax=100     |  S1: 15
    Csta=0, Cprop=0.     |  S2: 25
    DP => 9h30           |
         S1 : 10         |
         S1 : 10         |
                         |
    (pilotable) prod_1_1 |
    Pmin=20, Pmax=100    |
    Csta=100k, Cprop=100 |
    =#
    @testset "tso_capping_limitables_due_to_pilotable_pmin" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T09:30:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
                                                Dates.Second(90*60), Dates.Second(90*60))
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                20., 100.,
                                                100000., 100.,
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
                    SortedDict("wind_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.TO_DECIDE),
                               "prod_1_1" => SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.FREE))
                    )
        # initial generators state :
        generators_init_state = SortedDict("prod_1_1" => PSCOPF.OFF)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                        generators_init_state,
                                        uncertainties, nothing)
        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), Dates.DateTime("2015-01-01T07:00:00"), SortedDict(
            "wind_1_1" => PSCOPF.GeneratorSchedule("wind_1_1",
                SortedDict(),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>10.,"S2"=>10.)))
            ),
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{PSCOPF.GeneratorState}(PSCOPF.ON,
                                                                                        SortedDict("S1"=>PSCOPF.OFF, "S2"=>PSCOPF.ON))),
                SortedDict(Dates.DateTime("2015-01-01T11:00:00") => PSCOPF.UncertainValue{Float64}(missing,
                                                                                        SortedDict("S1"=>0.,"S2"=>20.)))
                ),
            )
        )

        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # indicates using slacks for feasibility (because of scenario S1)
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        # In S1 : Load=15, wind provides 10 => still missing 5 but pmin=20
        # => we want to reduce consumption by 5
        # But limit links the two scenarios and we can only produce to limit or uncertainty levels
        @test 5. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "prod_1_1", TS[1], "S1") < 1e-09
        @test 10. ≈ value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S1"])

        # In S2 : Load=25, wind provides 10 => still missing 15 but pmin=20
        # pilotable produces 20 => 5 extra prod (20+10 - 25) => reduce wind by 5.
        @test 5. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S2")
        @test 20. ≈ PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "prod_1_1", TS[1], "S2")
        @test -5. ≈ ( PSCOPF.get_prod_value(PSCOPF.get_tso_schedule(context), "wind_1_1", TS[1], "S2")
                    - PSCOPF.get_uncertainties(uncertainties[ech], "wind_1_1", TS[1], "S2") )
        @test value(result.lol_model.p_loss_of_load["bus_1", TS[1], "S2"]) < 1e-09
    end


    #=
    Capped power is the difference between the available power in a limitable
     and the decided production allowed for that limitable unit.

    TS: [11h]
    S: [S1]
                      bus 1
                        |
    (limitable) wind_1_1|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=2     |     load_1
      S1: 50            | S1: 15
                        |
    (limitable) wind_1_2|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |
      S1: 10            |
                        |
    (limitable) wind_1_3|
    Pmin=0, Pmax=100    |
    Csta=0, Cprop=1     |
      S1: 0             |
    =#
    @testset "tso_capping_in_schedule" begin
        TS = [DateTime("2015-01-01T11:00:00")]
        ech = DateTime("2015-01-01T07:00:00")
        network = PSCOPF.Networks.Network()
        PSCOPF.Networks.add_new_bus!(network, "bus_1")
        # Limitables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 2.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "wind_1_2", PSCOPF.Networks.LIMITABLE,
                                                0., 100.,
                                                0., 1.,
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
        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # Limitable injections
        @test 15. ≈ PSCOPF.get_prod_value(context.tso_schedule, "wind_1_1", TS[1], "S1")
        @test PSCOPF.get_prod_value(context.tso_schedule, "wind_1_2", TS[1], "S1") <= 1e-09
        @test PSCOPF.get_prod_value(context.tso_schedule, "wind_1_3", TS[1], "S1") <= 1e-09

        # Capped power
        @test (50. - 15) ≈ context.tso_schedule.capping["wind_1_1", TS[1], "S1"]
        @test (10. - 0)  ≈ context.tso_schedule.capping["wind_1_2", TS[1], "S1"]
        @test (0.  - 0)  ≈ context.tso_schedule.capping["wind_1_3", TS[1], "S1"]

        #Cost
        @test value(result.objective_model.start_cost) < 1e-09
        @test value(result.objective_model.prop_cost) ≈ (
              ((50. - 15.) * 2 + (10. - 0.) * 1. + (0. - 0.) * 1.)
        )
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
                        |        5             |
                        |                      |
    load                |                      |load
      S1: 180           |                      |  S1: 20

    available prod : 160
    total load : 200
    => need to cut 40
    optimization decides on how much load to shed on each bus:
        bus1: cuts at least 15 (cause branch is limited to 5 : 20-5 = 15)
        bus2: cuts at most 40-15=25
        bus1+bus2 = 40
    =#
    @testset "tso_cutting_load_by_bus_in_schedule" begin
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
        tso = PSCOPF.TSOOutFO()
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)
        PSCOPF.update_tso_schedule!(context, ech, result, firmness, tso)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_HAS_SLACK

        # possible production is 160 but load is 200
        @test 160. ≈ PSCOPF.get_prod_value(context.tso_schedule, "wind_1_1", TS[1], "S1")
        @test 40. ≈ (context.tso_schedule.loss_of_load_by_bus["bus_1", TS[1], "S1"]
                    + context.tso_schedule.loss_of_load_by_bus["bus_2", TS[1], "S1"])
        # cut conso by bus
        @test context.tso_schedule.loss_of_load_by_bus["bus_1", TS[1], "S1"] <= 25
        @test context.tso_schedule.loss_of_load_by_bus["bus_2", TS[1], "S1"] >= 15
    end

    #TODO : illustrate ptdf effect

end
