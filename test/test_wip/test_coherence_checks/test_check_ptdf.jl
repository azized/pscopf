using PSCOPF

using Test
using DataStructures

@testset "test_check_ptdf" begin

    network = PSCOPF.Networks.Network()
    PSCOPF.add_new_buses!(network, ["bus_1", "bus_2", "bus_3"])
    PSCOPF.add_new_branch!(network, "branch_1_2", 1000.)
    PSCOPF.add_new_branch!(network, "branch_1_3", 1000.)

    @testset "definition_example" begin
        ptdf = PSCOPF.Networks.PTDFValues(
                        "branch_1_2" => SortedDict("bus_1"=>0.,
                                                    "bus_2"=>0.,
                                                    "bus_3"=>0.),
                        "branch_1_3" => SortedDict("bus_1"=>0.,
                                                    "bus_2"=>0.,
                                                    "bus_3"=>0.)
        )

        @test PSCOPF.check_ptdf(ptdf, network)
    end

    @testset "missing_branch_1_3_entry" begin
        ptdf = SortedDict{String,SortedDict{String, Float64}}(
                        "branch_1_2" => SortedDict("bus_1"=>0.,
                                                    "bus_2"=>0.,
                                                    "bus_3"=>0.)
        )

        @test !PSCOPF.check_ptdf(ptdf, network)
    end

    @testset "missing_bus_2_in_branch_1_2" begin
        ptdf = SortedDict{String,SortedDict{String, Float64}}(
                        "branch_1_2" => SortedDict("bus_1"=>0.,
                                                    "bus_3"=>0.),
                        "branch_1_3" => SortedDict("bus_1"=>0.,
                                                    "bus_2"=>0.,
                                                    "bus_3"=>0.)
        )

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