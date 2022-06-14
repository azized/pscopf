using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_bilevel_pilotables_costing" begin

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")
    function create_instance(load_1, load_2,
                             prop_cost_1, prop_cost_2,
                            limit::Float64=35.,
                            logs=nothing)
        network = PSCOPFFixtures.network_2buses(limit=limit)
        # Pilotables
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "prod_1_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., prop_cost_1,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., prop_cost_2,
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
      S1: 10            |                      | S1: 100

    prod_1_1 is cheaper than prod_2_1
     but using only prod_1_1, like a market would, will cause RSO problems
    => The TSO needs to oblige
            prod_1_1 to produce at most 10+35=45 => P_1_1 in [0,45] => cost : 200-45 = 1*155
            prod_2_1 to produce at least 100-35=65 => P_2_1 in [65,200] => cost : 65 = 1*65
        Imposing one of the preceding conditions guarantees RSO constraints
        => TSO will impose one of the two preceding conditions, the cheaper one => option 2

    =#
    @testset "configured_same_cost" begin

        context = create_instance(10., 100.,
                                10., 50.,
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

        #Market chooses the levels sets the bounds for pilotable production
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 65. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test objective_value(result.upper.model) ≈ ( 65 )

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

    TSOBilevelConfigs.USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST allows using the units' prop cost as
     coeff in the objective function for limiting pilotables.

    prod_1_1 is cheaper than prod_2_1
     but using only prod_1_1, like a market would, will cause RSO problems
    => The TSO needs to oblige
            bus1 to produce at most 10+35=45 => P_1_1 in [0,45] => cost : 200-45 = 10*155 = 1550
            bus2 to produce at least 100-35=65 => P_2_1 in [65,200] => cost : 65 = 50*65 = 3250
        Imposing one of the preceding conditions guarantees RSO constraints
        => TSO will impose one of the two preceding conditions, the cheaper one => option 1

    =#
    @testset "use_prop_costs" begin

        context = create_instance(10., 100.,
                                10., 50.,
                                35.)

        tso = PSCOPF.TSOBilevel(PSCOPF.TSOBilevelConfigs(USE_UNITS_PROP_COST_AS_TSO_BOUNDING_COST=true))

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
        @test 45. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_1_1",TS[1],"S1"])
        #prod_2_1 bounding + EOD, guarantees RSO constraints
        @test value(result.upper.pilotable_model.p_imposition_min["prod_2_1",TS[1],"S1"]) < 1e-09
        @test 200. ≈ value(result.upper.pilotable_model.p_imposition_max["prod_2_1",TS[1],"S1"])

        #Market chooses the levels sets the bounds for pilotable production
        @test 45. ≈ value(result.lower.pilotable_model.p_injected["prod_1_1",TS[1],"S1"])
        @test 65. ≈ value(result.lower.pilotable_model.p_injected["prod_2_1",TS[1],"S1"])

        @test objective_value(result.upper.model) ≈ 1550.

    end

end
