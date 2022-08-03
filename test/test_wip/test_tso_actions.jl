using PSCOPF

using Test
using Dates

@testset verbose=true "test_tso_actions" begin

    ts = Dates.DateTime("2015-01-01T14:00:00")
    ts2 = Dates.DateTime("2015-01-01T14:30:00")

    @testset "empty" begin
        tso_actions = PSCOPF.TSOActions()

        @test isempty( PSCOPF.get_limitations(tso_actions) )
        @test ismissing( PSCOPF.get_limitation(tso_actions, "lim_1", ts) )

        @test isempty( PSCOPF.get_impositions(tso_actions) )
        @test ismissing( PSCOPF.get_imposition(tso_actions, "imp", ts, "S1") )
    end

    @testset "definitive_value" begin
        tso_actions = PSCOPF.TSOActions()
        PSCOPF.set_limitation_definitive_value!(tso_actions, "lim_1", ts, 300.)
        PSCOPF.set_limitation_value!(tso_actions, "lim_2", ts, "S1", 91.)
        PSCOPF.set_limitation_value!(tso_actions, "lim_2", ts, "S2", 92.)
        PSCOPF.set_imposition_value!(tso_actions, "imp", ts, "S1", 2., 10.)
        PSCOPF.set_imposition_value!(tso_actions, "imp", ts2, "S1", 55., 55.)

        @test length( PSCOPF.get_limitations(tso_actions) ) == 2
        @test PSCOPF.get_limitation(tso_actions, "lim_1", ts) == 300.
        @test PSCOPF.get_limitation(tso_actions, "lim_2", ts, "S1") == 91.
        @test PSCOPF.get_limitation(tso_actions, "lim_2", ts, "S2") == 92.

        @test length( PSCOPF.get_impositions(tso_actions) ) == 2
        @test PSCOPF.get_imposition(tso_actions, "imp", ts, "S1") == (2., 10.)
        @test PSCOPF.get_imposition(tso_actions, "imp", ts2, "S1") == (55., 55.)

        @test_throws ErrorException PSCOPF.get_imposition_level(tso_actions, "imp", ts, "S1") #cause imposition[1] !=  imposition[2] => bounds not a level
        @test PSCOPF.get_imposition_level(tso_actions, "imp", ts2, "S1") == 55.
    end

end
