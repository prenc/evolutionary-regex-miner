include("../src/evolution.jl")

using Test
using .evolution: remove_event, add_branch_or, crossover, find_brackets

@testset "Remove event mutation" begin
    @testset "abcde" begin
        @test "bcde" == remove_event("abcde", idx = 1)
        @test "acde" == remove_event("abcde", idx = 2)
        @test "abde" == remove_event("abcde", idx = 3)
        @test "abce" == remove_event("abcde", idx = 4)
        @test "abcd" == remove_event("abcde", idx = 5)
    end

    @testset "(a|b)" begin
        @test "b" == remove_event("(a|b)", idx = 2)
        @test "a" == remove_event("(a|b)", idx = 4)
    end

    @testset "(a+|b)+" begin
        @test "b" == remove_event("(a+|b)+", idx = 2)
        @test "a+" == remove_event("(a+|b)+", idx = 5)
    end

    @testset "(abc|b)+" begin
        @test "b" == remove_event("(abc|b)+", idx = 2)
        @test "(ac|b)+" == remove_event("(abc|b)+", idx = 3)
        @test "b" == remove_event("(abc|b)+", idx = 4)
        @test "abc" == remove_event("(abc|b)+", idx = 6)
    end

    @testset "[ab]{2}" begin
        @test "b" == remove_event("[ab]{2}", idx = 2)
        @test "a" == remove_event("[ab]{2}", idx = 3)
    end

    @testset "a[bcd]{2}e" begin
        @test "[bcd]{2}e" == remove_event("a[bcd]{2}e", idx = 1)
        @test "a[cd]{2}e" == remove_event("a[bcd]{2}e", idx = 3)
        @test "a[bd]{2}e" == remove_event("a[bcd]{2}e", idx = 4)
        @test "a[bc]{2}e" == remove_event("a[bcd]{2}e", idx = 5)
        @test "a[bcd]{2}" == remove_event("a[bcd]{2}e", idx = 10)
    end

    @testset "(a[bc]{2}|d)" begin
        @test "d" == remove_event("(a[bc]{2}|d)", idx = 2)
        @test "(ac|d)" == remove_event("(a[bc]{2}|d)", idx = 4)
        @test "(ab|d)" == remove_event("(a[bc]{2}|d)", idx = 5)
        @test "a[bc]{2}" == remove_event("(a[bc]{2}|d)", idx = 11)
    end

    @testset "(ab)+c" begin
        @test "b+c" == remove_event("(ab)+c", idx = 2)
        @test "a+c" == remove_event("(ab)+c", idx = 3)
        @test "(ab)+" == remove_event("(ab)+c", idx = 6)
    end
end

@testset "Add branch mutation" begin
    @testset "a" begin
        @test "(a|b)" == add_branch_or("a", "b", idx = 1)
    end

    @testset "a+" begin
        @test "(a|b)+" == add_branch_or("a+", "b", idx = 1)
    end

    @testset "a(b|c)d" begin
        @test "(a|e)(b|c)d" == add_branch_or("a(b|c)d", "e", idx = 1)
        @test "a((b|e)|c)d" == add_branch_or("a(b|c)d", "e", idx = 3)
        @test "a(b|(c|e))d" == add_branch_or("a(b|c)d", "e", idx = 5)
        @test "a(b|c)(d|e)" == add_branch_or("a(b|c)d", "e", idx = 7)
    end

    @testset "(ab|c)" begin
        @test "((a|e)b|c)" == add_branch_or("(ab|c)", "e", idx = 2)
        @test "(a(b|e)|c)" == add_branch_or("(ab|c)", "e", idx = 3)
        @test "(ab|(c|e))" == add_branch_or("(ab|c)", "e", idx = 5)
    end

    @testset "a[bc]{2}d" begin
        @test "(a|e)[bc]{2}d" == add_branch_or("a[bc]{2}d", "e", idx = 1)
        @test "a[bc]{2}(d|e)" == add_branch_or("a[bc]{2}d", "e", idx = 9)
    end
end

@testset "Crossover" begin
    @testset "a(b|c) and d(e|f)" begin
        @test "d(b|c)" == crossover("a(b|c)", "d(e|f)", idx1 = 1, idx2 = 1)
        @test "(e|f)(b|c)" == crossover("a(b|c)", "d(e|f)", idx1 = 1, idx2 = 2)
        @test "e(b|c)" == crossover("a(b|c)", "d(e|f)", idx1 = 1, idx2 = 3)
        @test "f(b|c)" == crossover("a(b|c)", "d(e|f)", idx1 = 1, idx2 = 5)

        @test "ad" == crossover("a(b|c)", "d(e|f)", idx1 = 2, idx2 = 1)

        @test "a((e|f)|c)" == crossover("a(b|c)", "d(e|f)", idx1 = 3, idx2 = 2)
    end

    @testset "a(b|c+)+ and d(e+|f)+" begin
        @test "d(b|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 1, idx2 = 1)
        @test "(e+|f)+(b|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 1, idx2 = 2)
        @test "e+(b|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 1, idx2 = 3)
        @test "f(b|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 1, idx2 = 6)
        @test "ad" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 2, idx2 = 1)
        @test "a(e+|f)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 2, idx2 = 2)
        @test "ae+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 2, idx2 = 3)
        @test "af" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 2, idx2 = 6)
        @test "a(d|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 3, idx2 = 1)
        @test "a((e+|f)+|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 3, idx2 = 2)
        @test "a(e+|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 3, idx2 = 3)
        @test "a(f|c+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 3, idx2 = 6)
        @test "a(b|d)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 5, idx2 = 1)
        @test "a(b|(e+|f)+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 5, idx2 = 2)
        @test "a(b|e+)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 5, idx2 = 3)
        @test "a(b|f)+" == crossover("a(b|c+)+", "d(e+|f)+", idx1 = 5, idx2 = 6)
    end

    @testset "(a|[bc]{2})d and (e+|f)+[gh]{2}" begin
        @test "(e+|f)+d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 1, idx2 = 1)
        @test "e+d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 1, idx2 = 2)
        @test "fd" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 1, idx2 = 5)
        @test "[gh]{2}d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 1, idx2 = 8)
        @test "((e+|f)+|[bc]{2})d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 2, idx2 = 1)
        @test "(e+|[bc]{2})d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 2, idx2 = 2)
        @test "(f|[bc]{2})d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 2, idx2 = 5)
        @test "([gh]{2}|[bc]{2})d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 2, idx2 = 8)
        @test "(a|(e+|f)+)d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 4, idx2 = 1)
        @test "(a|e+)d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 4, idx2 = 2)
        @test "(a|f)d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 4, idx2 = 5)
        @test "(a|[gh]{2})d" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 4, idx2 = 8)
        @test "(a|[bc]{2})(e+|f)+" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 12, idx2 = 1)
        @test "(a|[bc]{2})e+" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 12, idx2 = 2)
        @test "(a|[bc]{2})f" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 12, idx2 = 5)
        @test "(a|[bc]{2})[gh]{2}" == crossover("(a|[bc]{2})d", "(e+|f)+[gh]{2}", idx1 = 12, idx2 = 8)
    end
end

@testset "Find branch indicies function" begin
    @testset "a(b|c)" begin
        @test (nothing, nothing, nothing, nothing) == find_brackets("a(b|c)", 1)
        @test (2, 4, 6, false) == find_brackets("a(b|c)", 2)
        @test (2, 4, 6, false) == find_brackets("a(b|c)", 3)
        @test (2, 4, 6, false) == find_brackets("a(b|c)", 4)
        @test (2, 4, 6, false) == find_brackets("a(b|c)", 5)
        @test (2, 4, 6, false) == find_brackets("a(b|c)", 6)
    end

    @testset "a(b+|c)" begin
        @test (nothing, nothing, nothing, nothing) == find_brackets("a(b+|c)", 1)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 2)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 3)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 4)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 5)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 6)
        @test (2, 5, 7, false) == find_brackets("a(b+|c)", 7)
    end

    @testset "(a|(b|(c|d)))" begin
        @test (1, 3, 13, false) == find_brackets("(a|(b|(c|d)))", 1)
        @test (1, 3, 13, false) == find_brackets("(a|(b|(c|d)))", 2)
        @test (1, 3, 13, false) == find_brackets("(a|(b|(c|d)))", 3)
        @test (4, 6, 12, false) == find_brackets("(a|(b|(c|d)))", 4)
        @test (4, 6, 12, false) == find_brackets("(a|(b|(c|d)))", 5)
        @test (4, 6, 12, false) == find_brackets("(a|(b|(c|d)))", 6)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d)))", 7)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d)))", 8)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d)))", 9)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d)))", 10)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d)))", 11)
        @test (4, 6, 12, false) == find_brackets("(a|(b|(c|d)))", 12)
        @test (1, 3, 13, false) == find_brackets("(a|(b|(c|d)))", 13)
    end

    @testset "(a|(b|(c|d))+)" begin
        @test (1, 3, 14, false) == find_brackets("(a|(b|(c|d))+)", 1)
        @test (1, 3, 14, false) == find_brackets("(a|(b|(c|d))+)", 2)
        @test (1, 3, 14, false) == find_brackets("(a|(b|(c|d))+)", 3)
        @test (4, 6, 12, true) == find_brackets("(a|(b|(c|d))+)", 4)
        @test (4, 6, 12, true) == find_brackets("(a|(b|(c|d))+)", 5)
        @test (4, 6, 12, true) == find_brackets("(a|(b|(c|d))+)", 6)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d))+)", 7)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d))+)", 8)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d))+)", 9)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d))+)", 10)
        @test (7, 9, 11, false) == find_brackets("(a|(b|(c|d))+)", 11)
        @test (4, 6, 12, true) == find_brackets("(a|(b|(c|d))+)", 12)
        @test (4, 6, 12, true) == find_brackets("(a|(b|(c|d))+)", 13)
        @test (1, 3, 14, false) == find_brackets("(a|(b|(c|d))+)", 14)
    end

    @testset "((a|b)+|(c|[de]{2}))+" begin
        @test (1, 8, 20, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 1)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 2)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 3)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 4)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 5)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 6)
        @test (2, 4, 6, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 7)
        @test (1, 8, 20, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 8)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 9)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 10)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 11)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 12)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 13)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 14)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 15)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 16)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 17)
        @test (9, 11, 19, false) == find_brackets("((a|b)+|(c|[de]{2}))+", 19)
        @test (1, 8, 20, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 20)
        @test (1, 8, 20, true) == find_brackets("((a|b)+|(c|[de]{2}))+", 21)
    end

    @testset "(ab)+c" begin
        @test (1, nothing, 4, true) == find_brackets("(ab)+c", 1)
        @test (1, nothing, 4, true) == find_brackets("(ab)+c", 2)
        @test (1, nothing, 4, true) == find_brackets("(ab)+c", 3)
        @test (1, nothing, 4, true) == find_brackets("(ab)+c", 4)
        @test (1, nothing, 4, true) == find_brackets("(ab)+c", 5)
        @test (nothing, nothing, nothing, nothing) == find_brackets("(ab)+c", 6)
    end
end
