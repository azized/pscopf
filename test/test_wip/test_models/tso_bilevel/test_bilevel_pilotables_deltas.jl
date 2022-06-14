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
                                                0., 40.,
                                                0., 10.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_2", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., 12.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                0., 500.,
                                                0., 50.,
                                                Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_1", TS[1], "S1", load_1)
        PSCOPF.add_uncertainty!(uncertainties, ech, "bus_2", TS[1], "S1", load_2)

        gen_initial_state = SortedDict{String,PSCOPF.GeneratorState}(
                                "prod_1_1" => PSCOPF.ON,
                                "prod_1_2" => PSCOPF.OFF,
                                "prod_2_1" => PSCOPF.ON,
                            )
        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    gen_initial_state, uncertainties, nothing, logs)

        return context
    end


    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=40     |                      |Pmin=0, Pmax=500
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
    (limitable) prod_1_2|                      |
    Pmin=0, Pmax=200    |                      |
    Csta=0, Cprop=12    |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1:  10           |                      | S1: 500

      ideally the market would use prod_1_1 and prod_1_2 first
      but the RSO constraint would not be satisfied
      => TSO needs to impose a maximum prod of 45 on bus 1 or a minimum prod of 465 on bus 2
      cost of unit prod_2_1 is high => impositions are high too (FIXME? this behaviour needs changing wrt to maximum imposition)
      TSO will :
        impose a maximum prod of 45 on bus 1:
            force prod_1_1 to 0 (it's the cheapest to impose)
            maximum of 45 on prod_1_2
            free prod_2_1

      The market is then obliged to use prod_1_2 which is more expensive than prod_1_1,
       when both satisfy the RSO constraint

    =#
    @testset "no_reference_market_no_delta" begin

        context = create_instance(10., 500.,
                                35.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(CONSIDER_DELTAS=false, USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        # #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]) < 1e-09
        #prod_1_2 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_2",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_2",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 500. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_2",TS[1],"S1"])
        @test 465. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ ( 40. *10. + 155. *12. )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ ( 45*12 + 465*50. )

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=40     |                      |Pmin=0, Pmax=500
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
    (limitable) prod_1_2|                      |
    Pmin=0, Pmax=200    |                      |
    Csta=0, Cprop=12    |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1:  10           |                      | S1: 500

      refenece market schedule is not considered here in terms of deltas
        (it still affects the cost function if there were previous UC decisions).
        => same behaviour as "no_reference_market_no_delta"

      ideally the market would use prod_1_1 and prod_1_2 first
      but the RSO constraint would not be satisfied
      => TSO needs to impose a maximum prod of 45 on bus 1 or a minimum prod of 465 on bus 2
      cost of unit prod_2_1 is high => impositions are high too (FIXME? this behaviour needs changing wrt to maximum imposition)
      TSO will :
        impose a maximum prod of 45 on bus 1:
            force prod_1_1 to 0 (it's the cheapest to impose)
            maximum of 45 on prod_1_2
            free prod_2_1

      The market is then obliged to use prod_1_2 which is more expensive than prod_1_1,
       when both satisfy the RSO constraint
    =#
    @testset "reference_market_no_delta" begin

        context = create_instance(10., 500.,
                                35.)

        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>40.,)))
                    ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>20.,)))
                    ),
            "prod_2_1" => PSCOPF.GeneratorSchedule("prod_2_1",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>450.,)))
                    )
            )
        )

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(CONSIDER_DELTAS=false, USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        # #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"]) < 1e-09
        #prod_1_2 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_2",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_2",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 500. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_2",TS[1],"S1"])
        @test 465. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ ( 40. *10. + 155. *12. )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ ( 45*12 + 465*50. )

    end

    #=
    TS: [11h]
    S: [S1]
                        bus 1                   bus 2
                        |                      |
    (limitable) prod_1_1|       "1_2"          |prod_2_1
    Pmin=0, Pmax=40     |                      |Pmin=0, Pmax=500
    Csta=0, Cprop=10    |                      |Csta=0, Cprop=50
                        |----------------------|
                        |         35           |
                        |                      |
                        |                      |
    (limitable) prod_1_2|                      |
    Pmin=0, Pmax=200    |                      |
    Csta=0, Cprop=12    |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      S1:  10           |                      | S1: 500

      ideally the market would use prod_1_1 and prod_1_2 first

      an artificial (not optimal) schedule is given here :
        prod_1_1:40, prod_1_2:20, prod_2_1:450
      => production on bus 1 is problematic (causes RSO constraint)

      => TSO needs to impose a maximum prod of 45 on bus 1 or a minimum prod of 465 on bus 2
      The TSO tries to deviate the least from the arket refence schedule
      cost of unit prod_2_1 is high => impositions are high too (FIXME? this behaviour needs changing wrt to maximum imposition)
      TSO will :
        impose a maximum prod of 45 on bus 1:
            constrain prod_1_1 (it's the cheapest to impose) => max at 20 (=> delta of 15)
            maxi at 20 on prod_1_2 (delta of 0)
            free prod_2_1 (production will compensate the reduced 15MW)
    =#
    @testset "reference_market_and_consider_delta" begin

        context = create_instance(10., 500.,
                                35.)

        context.market_schedule = PSCOPF.Schedule(PSCOPF.Market(), ech, SortedDict(
            "prod_1_1" => PSCOPF.GeneratorSchedule("prod_1_1",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>40.,)))
                    ),
            "prod_1_2" => PSCOPF.GeneratorSchedule("prod_1_2",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>20.,)))
                    ),
            "prod_2_1" => PSCOPF.GeneratorSchedule("prod_2_1",
                    SortedDict(),
                    SortedDict(TS[1] => PSCOPF.UncertainValue{PSCOPF.Float64}(missing,
                                                                                        SortedDict("S1"=>450.,)))
                    )
            )
        )

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(CONSIDER_DELTAS=true, USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))
        firmness = PSCOPF.compute_firmness(tso,
                                            ech, next_ech,
                                            TS, context)
        result = PSCOPF.run(tso, ech, firmness,
                    PSCOPF.get_target_timepoints(context),
                    context)

        # Solution is optimal
        @test PSCOPF.get_status(result) == PSCOPF.pscopf_OPTIMAL

        # #TSO RSO constraints are OK
        @test value(result.upper.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.upper.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        # #Market EOD constraints are OK
        @test value(result.lower.limitable_model.p_global_capping[TS[1],"S1"]) < 1e-09
        @test value(result.lower.lol_model.p_global_loss_of_load[TS[1],"S1"]) < 1e-09

        #TSO sets the bounds for pilotable production respecting units' pmin and pmax
        #prod_1_1 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_1",TS[1],"S1"]) < 1e-09
        @test 25. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_1_2 is not bound
        @test value(result.upper.pilotable_model.p_imposition_min["prod_1_2",TS[1],"S1"]) < 1e-09
        @test 20. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_2",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 500. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels respecting the bounds for pilotable production
        @test 25. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 20. ≈ value(result.lower.pilotable_model.p_injected["prod_1_2",TS[1],"S1"])
        @test 465. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test value(PSCOPF.get_upper_obj_expr(result)) ≈ ( 15. *10. + 180. *12. )
        @test value(PSCOPF.get_lower_obj_expr(result)) ≈ ( 25. *10 + 20. *12. + 465*50. )

    end

end
