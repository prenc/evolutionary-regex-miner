include("../src/scoring.jl")

using Test
using .scoring: get_bracket_levels

@testset "Evaluate bracket levels function" begin
    @test [1] == get_bracket_levels("(a|b)")
    @test [1, 2] == get_bracket_levels("a(a|(da|b))c")
    @test [1, 2, 3] == get_bracket_levels("a(a|((d|c)a|b))c")
    @test [1, 1, 2, 2, 3] == get_bracket_levels("(a)+((d|a)a|((d|c)a|b))c")
    @test [1, 2, 3, 4, 5] == get_bracket_levels("(a(b|(c|(d|(e|f)))))+")
end
