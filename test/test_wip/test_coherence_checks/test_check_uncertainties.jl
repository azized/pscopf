using PSCOPF

using Test
using Dates
using DataStructures

@testset "test_check_TS" begin

    network = PSCOPF.Networks.Network()
    PSCOPF.add_new_bus!(network, "bus_1")
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "limitable_1", PSCOPF.Networks.LIMITABLE,
                                                0., 100., #pmin, pmax : Not concerned ? min is always 0, max is the limitation
                                                0., 10., #start_cost, prop_cost : start cost is always 0 ?
                                                Dates.Second(0), Dates.Second(0)) #dmo, dp : always 0. ?
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "pilotable_1", PSCOPF.Networks.PILOTABLE,
                                                10., 200., #pmin, pmax
                                                450., 30., #start_cost, prop_cost
                                                Dates.Second(3600), Dates.Second(0)) #dmo, dp

    TS = [DateTime("2015-01-01T14:00:00")]
    ECH = [DateTime("2015-01-01T07:00:00"), DateTime("2015-01-01T11:00:00")]

    @testset "definition_example" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                "S2" => 20.
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 0.,
                                                                "S2" => 16.
                                                                )
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                "S2" => 20.
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 10.,
                                                                "S2" => 16.
                                                                )
                                ) ,
                            ),
                        )

        @test PSCOPF.check_uncertainties(uncertainties, network)
    end

    @testset "all_entries_need_to_have_the_same_scenarios" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                "S2" => 20.
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 10.,
                                                                # "S2" => 16. #Missing Scenario
                                                                )
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                "S2" => 20.
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 10.,
                                                                "S2" => 16.
                                                                )
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_same_scenarios(uncertainties)
    end

    @testset "all_entries_need_to_have_the_same_timesteps" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                #different/missing timestep
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_same_timesteps(uncertainties)
    end

    # in the execution we allow having more TS than what we need !
    @testset "all_specified_ts_are_listed" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )
        @test PSCOPF.check_uncertainties_contain_ts(uncertainties, [DateTime("2015-01-01T14:00:00")])

        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                #DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                #DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                #DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ),
                            "limitable_1" =>SortedDict(
                                #DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.),
                                DateTime("2015-01-01T14:15:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )
        @test !PSCOPF.check_uncertainties_contain_ts(uncertainties,
                                                    [DateTime("2015-01-01T14:00:00"),DateTime("2015-01-01T14:15:00")])
    end

    @testset "all_values_must_be_>=0" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => -10., # negative value
                                                                )
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 15.,
                                                                )
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict(
                                                                "S1" => 10.,
                                                                )
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_values(uncertainties, network)
    end

    @testset "limitable_injection_uncertainties<=pmax" begin
        # Allow extra ECHs
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 115.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 110.) # exceeds pmax (100)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 115.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 90.)
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_values(uncertainties, network)
    end

    @testset "all_listed_ids_are_either_a_bus_or_a_limitable_of_the_network" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "bus_2" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "bus_2" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_values(uncertainties, network)
    end

    @testset "all_limitable_generators_are_listed_for_all_ech" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            # missing :
                            # "limitable_1" =>SortedDict(
                            #     DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                            #     ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )

        @test !PSCOPF.check_uncertainties_limitables(uncertainties, network)
    end

    @testset "all_buses_are_listed_for_all_ech" begin
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            # missing :
                            # "bus_1" => SortedDict(
                            #     DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                            #     ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )
        @test !PSCOPF.check_uncertainties_buses(uncertainties, network)
    end

    @testset "at_least_all_ech_are_described" begin
        # Allow extra ECHs
        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        # Extra ECH :
                        DateTime("2015-01-01T10:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        DateTime("2015-01-01T11:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        )

        @test PSCOPF.check_uncertainties_contains_ech(uncertainties, ECH)

        uncertainties = PSCOPF.Uncertainties(
                        DateTime("2015-01-01T07:00:00") => SortedDict(
                            "bus_1" => SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                                ),
                            "limitable_1" =>SortedDict(
                                DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                                ) ,
                            ),
                        # missing :
                        # DateTime("2015-01-01T11:00:00") => SortedDict(
                        #     "bus_1" => SortedDict(
                        #         DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 15.)
                        #         ),
                        #     "limitable_1" =>SortedDict(
                        #         DateTime("2015-01-01T14:00:00") => SortedDict("S1" => 10.)
                        #         ) ,
                        #     ),
                        )

        @test !PSCOPF.check_uncertainties_contains_ech(uncertainties, ECH)
    end

end