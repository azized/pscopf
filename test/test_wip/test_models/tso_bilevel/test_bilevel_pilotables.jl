using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_pilotables" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    function create_instance(load_1, load_2,
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., 10.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., 50.,
                                                Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", load_1)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", load_2)

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, nothing, logs)

        return context
    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=200    |                      |Pmin=0, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 30            |                      | S1: 30


    Market :
        The market uses prod_1_1 to satisfy all the demand
        This causes no RSO constraints

    TSO : does not to take any action
        lol_min = p_global_loss_of_load = 0.
        prod_1_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)
        prod_2_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)

    =#
    @testset "no_risk_of_breaking_rso_constraint" begin

        context = create_instance(30., 30.,
                                35.)

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
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test 60. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"]) < 1e-09

        @test value(PSCOPF.get_upper_obj_expr(result)) < 1e-09
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ (60. * 10.)

    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=200    |                      |Pmin=0, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 10            |                      | S1: 100

    Market :
        market would use prod_1_1 at 110MW cause it is cheaper
        This would cause a RSO constraint (flow of 100 > 35)

    TSO needs to anticipate to prevent that:
        option 1 : oblige prod_1_1 to produce at most 10+35=45 => P_1_1 in [0,45]
            => cost : 200-45=155
        option 2 : oblige prod_2_1 to produce at least 100-35=65 => P_2_1 in [65,200]
            => cost : 65

    Only obliging one of the options + EOD constraint is enough.
    The cheaper option is chosen : option 2
        lol_min = p_global_loss_of_load = 0.
        prod_1_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)
        prod_2_1 : p_imposition_min=65.  p_imposition_max=200

    Market solution :
        lol = p_loss_of_load = 0.
        prod_1_1 : 45.
        prod_2_1 : 65.
    =#
    @testset "rso_constraint_requires_setting_pilotables_bounds" begin

        context = create_instance(10., 100.,
                                35.)

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
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test 65. ≈ value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"])
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 65. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ ( 65 )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ ( 45*10 + 65*50. )

    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=200    |                      |Pmin=0, Pmax=200
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         500          |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 400           |                      | S1: 20

    Market :
        Possible prod : 400,  Demand : 420
        market will use prod_1_1 and prod_2_1 (=> 400MW)
        market will need to cut 20MW
        This will not cause RSO constraints violation

    TSO does not need to intervene
        lol_min = p_global_loss_of_load = 0.
        prod_1_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)
        prod_2_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)

        The TSO will locate the LoL on either of bus 1 or 2
        p_loss_of_load[bus1] + p_loss_of_load[bus2] = 20

    Market solution :
        lol = p_loss_of_load = 20.
        prod_1_1 : 200.
        prod_2_1 : 200.

    =#
    @testset "eod_problem_requires_cutting_conso" begin

        context = create_instance(400., 20.,
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

        #TSO does not need to bound pilotables
        #prod_1_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (20 *tso.configs.TSO_LOL_PENALTY)
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ ( 200*10. + 200*50.
                                                        + 20 *tso.configs.MARKET_LOL_PENALTY)

    end

    #=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|                      |prod_2_1
    Pmin=0, Pmax=200    |                      |Pmin=0, Pmax=200
    Csta=0, Cprop=10    |        "1_2"         |Csta=0, Cprop=50
                        |----------------------|
                        |         180          |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1: 400           |                      | S1: 20

    Market :
        Possible prod : 400,  Demand : 420
        market will use prod_1_1 and prod_2_1 (=> 400MW)
        market will need to cut 20MW

    If load is cut from bus 1 => no problem
    If load is cut from bus 2 => RSO constraint violation

    TSO needs to prevent this by locating the cut conso on bus 2
        option 1 : limit prod2 to 0-20
        lol_min = p_global_loss_of_load = 0.
        prod_1_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)
        prod_2_1 : p_imposition_min=0.  p_imposition_max=200 (unchanged bounds)

        The TSO will locate the LoL on bus 1 to avoid the constraint violation
        p_loss_of_load[bus1] = 20
        p_loss_of_load[bus2] = 0

    Market solution :
        lol = p_loss_of_load = 20.
        prod_1_1 : 200.
        prod_2_1 : 200.

    =#
    @testset "eod_and_rso_problems_require_cutting_conso_and_locating_it" begin

        context = create_instance(400., 20.,
                                180.)

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

        # RSO locates cut conso : due to flow limit, conso must be cut on bus_1
        @test 20. ≈ value(result.upper.lol_model.p_loss_of_load["bus_1", TS[1],"S1"])
        @test value(result.upper.lol_model.p_loss_of_load["bus_2", TS[1],"S1"]) < 1e-09

        #TSO does not need to bound pilotables
        #prod_1_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 200. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ (20 * tso.configs.TSO_LOL_PENALTY) #LoL due to Market

    end

end
