include("../src/evolution.jl")

using evolution: find_brackets

@testset "Find brackets function" begin
    @testset "a(b|c)" begin
        regex = "a(b|c)"
        @test (nothing, nothing, nothing, nothing) == find_brackets(regex, 1)
        @test (2, 4, 6, false) == find_brackets(regex, 2)
        @test (2, 4, 6, false) == find_brackets(regex, 3)
        @test (2, 4, 6, false) == find_brackets(regex, 4)
        @test (2, 4, 6, false) == find_brackets(regex, 5)
        @test (2, 4, 6, false) == find_brackets(regex, 6)
    end

    @testset "a(b+|c)" begin
        regex = "a(b+|c)"
        @test (nothing, nothing, nothing, nothing) == find_brackets(regex, 1)
        @test (2, 5, 7, false) == find_brackets(regex, 2)
        @test (2, 5, 7, false) == find_brackets(regex, 3)
        @test (2, 5, 7, false) == find_brackets(regex, 4)
        @test (2, 5, 7, false) == find_brackets(regex, 5)
        @test (2, 5, 7, false) == find_brackets(regex, 6)
        @test (2, 5, 7, false) == find_brackets(regex, 7)
    end

    @testset "(a|(b|(c|d)))" begin
        regex = "(a|(b|(c|d)))"
        @test (1, 3, 13, false) == find_brackets(regex, 1)
        @test (1, 3, 13, false) == find_brackets(regex, 2)
        @test (1, 3, 13, false) == find_brackets(regex, 3)
        @test (4, 6, 12, false) == find_brackets(regex, 4)
        @test (4, 6, 12, false) == find_brackets(regex, 5)
        @test (4, 6, 12, false) == find_brackets(regex, 6)
        @test (7, 9, 11, false) == find_brackets(regex, 7)
        @test (7, 9, 11, false) == find_brackets(regex, 8)
        @test (7, 9, 11, false) == find_brackets(regex, 9)
        @test (7, 9, 11, false) == find_brackets(regex, 10)
        @test (7, 9, 11, false) == find_brackets(regex, 11)
        @test (4, 6, 12, false) == find_brackets(regex, 12)
        @test (1, 3, 13, false) == find_brackets(regex, 13)
    end

    @testset "(a|(b|(c|d))+)" begin
        regex = "(a|(b|(c|d))+)"
        @test (1, 3, 14, false) == find_brackets(regex, 1)
        @test (1, 3, 14, false) == find_brackets(regex, 2)
        @test (1, 3, 14, false) == find_brackets(regex, 3)
        @test (4, 6, 12, true) == find_brackets(regex, 4)
        @test (4, 6, 12, true) == find_brackets(regex, 5)
        @test (4, 6, 12, true) == find_brackets(regex, 6)
        @test (7, 9, 11, false) == find_brackets(regex, 7)
        @test (7, 9, 11, false) == find_brackets(regex, 8)
        @test (7, 9, 11, false) == find_brackets(regex, 9)
        @test (7, 9, 11, false) == find_brackets(regex, 10)
        @test (7, 9, 11, false) == find_brackets(regex, 11)
        @test (4, 6, 12, true) == find_brackets(regex, 12)
        @test (4, 6, 12, true) == find_brackets(regex, 13)
        @test (1, 3, 14, false) == find_brackets(regex, 14)
    end

    @testset "((a|b)+|(c|d))+" begin
        regex = "((a|b)+|(c|d))+"
        @test (1, 8, 14, true) == find_brackets(regex, 1)
        @test (2, 4, 6, true) == find_brackets(regex, 2)
        @test (2, 4, 6, true) == find_brackets(regex, 3)
        @test (2, 4, 6, true) == find_brackets(regex, 4)
        @test (2, 4, 6, true) == find_brackets(regex, 5)
        @test (2, 4, 6, true) == find_brackets(regex, 6)
        @test (2, 4, 6, true) == find_brackets(regex, 7)
        @test (1, 8, 14, true) == find_brackets(regex, 8)
        @test (9, 11, 13, false) == find_brackets(regex, 9)
        @test (9, 11, 13, false) == find_brackets(regex, 10)
        @test (9, 11, 13, false) == find_brackets(regex, 11)
        @test (9, 11, 13, false) == find_brackets(regex, 12)
        @test (9, 11, 13, false) == find_brackets(regex, 13)
        @test (1, 8, 14, true) == find_brackets(regex, 14)
        @test (1, 8, 14, true) == find_brackets(regex, 15)
    end
end
