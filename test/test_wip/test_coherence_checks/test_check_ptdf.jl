using PSCOPF

using Test
using DataStructures

@testset "test_check_ptdf" begin

    network = PSCOPF.Networks.Network()
    PSCOPF.add_new_buses!(network, ["bus_1", "bus_2", "bus_3"])
    PSCOPF.add_new_branch!(network, "branch_1_2", 1000.)
    PSCOPF.add_new_branch!(network, "branch_1_3", 1000.)

    @testset "ptdf_dict" begin
        @testset "working_ptdf" begin
            ptdf = PSCOPF.Networks.PTDFDict(
                "BASECASE" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_2" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_3" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
            )

            @test PSCOPF.check_ptdf_case_entries(ptdf, network)
            @test PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "ptdf_with_missing_n-1" begin
            ptdf = PSCOPF.Networks.PTDFDict(
                "BASECASE" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_2" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
            )

            @test_broken !PSCOPF.check_ptdf_case_entries(ptdf, network)
            @test_broken !PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "ptdf_with_missing_branch_values" begin
            ptdf = PSCOPF.Networks.PTDFDict(
                "BASECASE" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    # "branch_1_3" => SortedDict("bus_1"=>0.,
                                    #                             "bus_2"=>0.,
                                    #                             "bus_3"=>0.)
                ),
                "branch_1_2" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_3" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
            )

            @test PSCOPF.check_ptdf_case_entries(ptdf, network)
            @test !PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "ptdf_with_missing_bus_value" begin
            ptdf = PSCOPF.Networks.PTDFDict(
                "BASECASE" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_2" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
                "branch_1_3" => PSCOPF.Networks.PTDFValues(
                                    "branch_1_2" => SortedDict("bus_1"=>0.,
                                                                "bus_2"=>0.,
                                                                "bus_3"=>0.),
                                    "branch_1_3" => SortedDict("bus_1"=>0.,
                                                                # "bus_2"=>0.,
                                                                "bus_3"=>0.)
                ),
            )

            @test PSCOPF.check_ptdf_case_entries(ptdf, network)
            @test !PSCOPF.check_ptdf(ptdf, network)
        end

    end

    @testset "check_ptdf_values" begin
        @testset "definition_example" begin
            ptdf = PSCOPF.Networks.PTDFValues(
                            "branch_1_2" => SortedDict("bus_1"=>0.,
                                                        "bus_2"=>0.,
                                                        "bus_3"=>0.),
                            "branch_1_3" => SortedDict("bus_1"=>0.,
                                                        "bus_2"=>0.,
                                                        "bus_3"=>0.)
            )

            @test PSCOPF.check_ptdf_branch_entries(ptdf, network)
            @test PSCOPF.check_ptdf_bus_entries(ptdf, network)
            @test PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "missing_branch_1_3_entry" begin
            ptdf = PSCOPF.Networks.PTDFValues(
                            "branch_1_2" => SortedDict("bus_1"=>0.,
                                                        "bus_2"=>0.,
                                                        "bus_3"=>0.)
            )

            @test !PSCOPF.check_ptdf_branch_entries(ptdf, network)
            @test PSCOPF.check_ptdf_bus_entries(ptdf, network)
            @test !PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "missing_bus_2_in_branch_1_2" begin
            ptdf = PSCOPF.Networks.PTDFValues(
                            "branch_1_2" => SortedDict("bus_1"=>0.,
                                                        "bus_3"=>0.),
                            "branch_1_3" => SortedDict("bus_1"=>0.,
                                                        "bus_2"=>0.,
                                                        "bus_3"=>0.)
            )

            @test PSCOPF.check_ptdf_branch_entries(ptdf, network)
            @test !PSCOPF.check_ptdf_bus_entries(ptdf, network)
            @test !PSCOPF.check_ptdf(ptdf, network)
        end

        @testset "missing_value_for_branch_1_2_bus_3" begin
            ptdf = SortedDict(
                                "branch_1_2" => SortedDict("bus_1"=>0.,
                                                            "bus_2"=>0.,
                                                            "bus_3"=>missing),
                                "branch_1_3" => SortedDict("bus_1"=>0.,
                                                            "bus_2"=>0.,
                                                            "bus_3"=>0.)
                            )
            @test_throws MethodError PSCOPF.check_ptdf(ptdf, network)
        end
    end

end