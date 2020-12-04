include("../src/scoring.jl")

using Test
using .scoring: get_bracket_levels

@testset "Get bracket levels function" begin
    regex = "(a|b)"
    @test [1] == get_bracket_levels(regex)
    regex = "a(a|(da|b))c"
    @test [1, 2] == get_bracket_levels(regex)
    regex = "a(a|((d|c)a|b))c"
    @test [1, 2, 3] == get_bracket_levels(regex)
    regex = "(a)+((d|a)a|((d|c)a|b))c"
    @test [1, 1, 2, 2, 3] == get_bracket_levels(regex)
    regex = "(a(b|(c|(d|(e|f)))))+"
    @test [1, 2, 3, 4, 5] == get_bracket_levels(regex)
end
