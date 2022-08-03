using PSCOPF

using Test
using Dates
using DataStructures

@testset "test_check_initial_state" begin

    network = PSCOPF.Networks.Network()
    PSCOPF.add_new_bus!(network, "bus_1")
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "limitable_1", PSCOPF.Networks.LIMITABLE,
                                                0., 10., #pmin, pmax : Not concerned ? min is always 0, max is the limitation
                                                0., 10., #start_cost, prop_cost : start cost is always 0 ?
                                                Dates.Second(0), Dates.Second(0)) #dmo, dp : always 0. ?
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "pilotable_1", PSCOPF.Networks.PILOTABLE,
                                                10., 200., #pmin, pmax
                                                450., 20., #start_cost, prop_cost
                                                Dates.Second(3600), Dates.Second(0)) #dmo, dp
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "pilotable_2", PSCOPF.Networks.PILOTABLE,
                                                20., 200., #pmin, pmax
                                                450., 20., #start_cost, prop_cost
                                                Dates.Second(3600), Dates.Second(0)) #dmo, dp
    PSCOPF.Networks.add_new_generator_to_bus!(network, "bus_1", "pilotable_3", PSCOPF.Networks.PILOTABLE,
                                                0., 200., #pmin, pmax
                                                0., 30., #start_cost, prop_cost
                                                Dates.Second(3600), Dates.Second(0)) #dmo, dp

    @testset "definition_example" begin
        initial_state = SortedDict(
                            "pilotable_1" => PSCOPF.ON,
                            "pilotable_2" => PSCOPF.ON,
        )
        @test PSCOPF.check_initial_state(initial_state, network)

        initial_state = SortedDict(
                            "pilotable_1" => PSCOPF.OFF,
                            "pilotable_2" => PSCOPF.ON,
                            "pilotable_3" => PSCOPF.ON,
                            "limitable_1" => PSCOPF.ON,
        )
        @test PSCOPF.check_initial_state(initial_state, network)
    end

    @testset "gen_with_pmin=0_can_be_ommited" begin
        # generators with pmin = 0 can be ommited (=> the case of limitables as well)
        initial_state = SortedDict(
                            "pilotable_1" => PSCOPF.OFF,
                            "pilotable_2" => PSCOPF.ON,
                            # "pilotable_3" => PSCOPF.ON,
                            # "limitable_1" => PSCOPF.ON,
        )
        @test PSCOPF.check_initial_state(initial_state, network)
    end

    @testset "gen_with_pmin=0_must_be_on" begin
        initial_state = SortedDict(
                            "pilotable_1" => PSCOPF.OFF,
                            "pilotable_2" => PSCOPF.ON,
                            "pilotable_3" => PSCOPF.OFF,
                            # "limitable_1" => PSCOPF.ON,
        )
        @test !PSCOPF.check_initial_state(initial_state, network)
    end

    @testset "all_gen_with_pmin>0_need_to_be_listed" begin
        initial_state = SortedDict(
                            "pilotable_1" => PSCOPF.OFF,
                            # "pilotable_2" => PSCOPF.ON,
                            "pilotable_3" => PSCOPF.ON,
                            "limitable_1" => PSCOPF.ON,
        )
        @test !PSCOPF.check_initial_state(initial_state, network)
    end

end