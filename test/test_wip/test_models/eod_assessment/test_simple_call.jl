using PSCOPF

using .PSCOPFFixtures

using Test
using JuMP
using Dates
using DataStructures
using Printf

@testset verbose=true "test_eod_assessment_call" begin

    function create_instance(limit_1, impositions_1, impositions_2,
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
                                                0., 200.,
                                                0., 10.,
                                                Dates.Second(0), Dates.Second(0))
        PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_2", "prod_2_1", PSCOPF.Networks.PILOTABLE,
                                                0., 200.,
                                                0., 50.,
                                                Dates.Second(0), Dates.Second(0))

        # Uncertainties
        uncertainties = PSCOPF.Uncertainties()

        #Fixme to args
        assessment_uncertainties = SortedDict("bus_1" => (30, 50),
                                            "bus_2" => (20, 70),
                                            "wind_1_1" => (50, 90)
                                            )

        context = PSCOPF.PSCOPFContext(network, TS, PSCOPF.PSCOPF_MODE_1,
                                    SortedDict{String,PSCOPF.GeneratorState}(), #gen_initial_state
                                    uncertainties, assessment_uncertainties, logs)

        PSCOPF.set_imposition_definitive_value!(PSCOPF.get_tso_actions(context), "prod_1_1", TS[1], impositions_1[1], impositions_1[2])
        PSCOPF.set_imposition_definitive_value!(PSCOPF.get_tso_actions(context), "prod_2_1", TS[1], impositions_2[1], impositions_2[2])
        PSCOPF.set_limitation_definitive_value!(PSCOPF.get_tso_actions(context), "wind_1_1", TS[1], limit_1)


        return context
    end

    TS = [DateTime("2015-01-01T11:00:00")]
    ech = DateTime("2015-01-01T07:00:00")
    next_ech = DateTime("2015-01-01T07:30:00")

    #=
                    bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=200    |                      |
        Csta=0, Cprop=1     |                      |
        [50-90]             |----------------------|
                            |         35           |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=20, Pmax=200   |                      | Pmin=20, Pmax=200
        Csta=10k, Cprop=10  |                      | Csta=50k, Cprop=50
        [0-0]               |                      | [20-200]
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        [30-50]             |                      | [20-70]
    =#
    @testset verbose=true "test_no_violation" begin
        limit_1 = 75.
        impositions_1 = (0., 0.)
        impositions_2 = (20., 200.)

        context = create_instance(limit_1, impositions_1, impositions_2,
                            35.,)

        assessment = PSCOPF.EODAssessment()
        result = PSCOPF.run(assessment, ech, TS, context)

        @test PSCOPF.is_validated(result)
    end

    #=
                    bus 1                   bus 2
                            |                      |
        (limitable) wind_1_1|       "1_2"          |
        Pmin=0, Pmax=200    |                      |
        Csta=0, Cprop=1     |                      |
        [50-90]             |----------------------|
                            |         35           |
                            |                      |
        (pilotable) prod_1_1|                      |(pilotable) prod_2_1
        Pmin=20, Pmax=200   |                      | Pmin=20, Pmax=200
        Csta=10k, Cprop=10  |                      | Csta=50k, Cprop=50
        [0-0]               |                      | [20-60]
                            |                      |
               load(bus_1)  |                      |load(bus_2)
        [30-50]             |                      | [20-70]
    =#
    @testset verbose=true "test_violation" begin
        limit_1 = 75.
        impositions_1 = (0., 0.)
        impositions_2 = (20., 60.)

        context = create_instance(limit_1, impositions_1, impositions_2,
                            35.,)

        assessment = PSCOPF.EODAssessment()
        result = PSCOPF.run(assessment, ech, TS, context)

        @test !PSCOPF.is_validated(result)
    end

end