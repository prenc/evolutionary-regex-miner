include("../src/evolution.jl")

using Test
using .evolution: add_event,
        remove_event,
        add_branch_or,
        remove_branch_or,
        add_branch_and,
        remove_branch_and,
        add_loop,
        remove_loop,
        pull_out,
        crossover,
        find_brackets

@testset "Add event mutation" begin
    @testset "Elementary addition" begin
        @test "ab" == add_event("b", "a", idx = 0)
        @test "ba" == add_event("b", "a", idx = 1)
    end

    @testset "Addition inside and branch" begin
        @test "[abc]{2}" == add_event("[ab]{2}", "c", idx = 1)
        @test "[abc]{2}" == add_event("[ab]{2}", "c", idx = 2)
        @test "[abc]{2}" == add_event("[ab]{2}", "c", idx = 3)
    end

    @testset "Addition inside or branch" begin
        @test "d(a|(b|c))" == add_event("(a|(b|c))", "d", idx = 0)
        @test "(da|(b|c))" == add_event("(a|(b|c))", "d", idx = 1)
        @test "(ad|(b|c))" == add_event("(a|(b|c))", "d", idx = 2)
        @test "(a|d(b|c))" == add_event("(a|(b|c))", "d", idx = 3)
        @test "(a|(db|c))" == add_event("(a|(b|c))", "d", idx = 4)
        @test "(a|(bd|c))" == add_event("(a|(b|c))", "d", idx = 5)
        @test "(a|(b|dc))" == add_event("(a|(b|c))", "d", idx = 6)
        @test "(a|(b|cd))" == add_event("(a|(b|c))", "d", idx = 7)
        @test "(a|(b|c)d)" == add_event("(a|(b|c))", "d", idx = 8)
        @test "(a|(b|c))d" == add_event("(a|(b|c))", "d", idx = 9)
    end
end

@testset "Remove event mutation" begin
    @testset "Elementary deletion" begin
        @test "bcde" == remove_event("abcde", idx = 1)
        @test "acde" == remove_event("abcde", idx = 2)
        @test "abde" == remove_event("abcde", idx = 3)
        @test "abce" == remove_event("abcde", idx = 4)
        @test "abcd" == remove_event("abcde", idx = 5)
    end

    @testset "Deletion in or branch" begin
        @test "b" == remove_event("(a|b)", idx = 2)
        @test "a" == remove_event("(a|b)", idx = 4)
        @test "b" == remove_event("(a+|b)+", idx = 2)
        @test "a+" == remove_event("(a+|b)+", idx = 5)

        @test "b" == remove_event("(abc|b)+", idx = 2)
        @test "(ac|b)+" == remove_event("(abc|b)+", idx = 3)
        @test "b" == remove_event("(abc|b)+", idx = 4)
        @test "abc" == remove_event("(abc|b)+", idx = 6)
    end

    @testset "Deletion inside and branch" begin
        @test "b" == remove_event("[ab]{2}", idx = 2)
        @test "a" == remove_event("[ab]{2}", idx = 3)

        @test "[bcd]{2}e" == remove_event("a[bcd]{2}e", idx = 1)
        @test "a[cd]{2}e" == remove_event("a[bcd]{2}e", idx = 3)
        @test "a[bd]{2}e" == remove_event("a[bcd]{2}e", idx = 4)
        @test "a[bc]{2}e" == remove_event("a[bcd]{2}e", idx = 5)
        @test "a[bcd]{2}" == remove_event("a[bcd]{2}e", idx = 10)
    end

    @testset "Deletion of and branch nested in or branch" begin
        @test "d" == remove_event("(a[bc]{2}|d)", idx = 2)
        @test "(ac|d)" == remove_event("(a[bc]{2}|d)", idx = 4)
        @test "(ab|d)" == remove_event("(a[bc]{2}|d)", idx = 5)
        @test "a[bc]{2}" == remove_event("(a[bc]{2}|d)", idx = 11)
    end

    @testset "Deletion of loops" begin
        @test "b+c" == remove_event("(ab)+c", idx = 2)
        @test "a+c" == remove_event("(ab)+c", idx = 3)
        @test "(ab)+" == remove_event("(ab)+c", idx = 6)
    end
end

@testset "Add or branch mutation" begin
    @testset "Elementary or branch addition" begin
        @test "(a|b)" == add_branch_or("a", "b", idx = 1)
        @test "(a|b)+" == add_branch_or("a+", "b", idx = 1)
    end

    @testset "Addition if nested or branches" begin
        @test "(a|e)(b|c)d" == add_branch_or("a(b|c)d", "e", idx = 1)
        @test "a((b|e)|c)d" == add_branch_or("a(b|c)d", "e", idx = 3)
        @test "a(b|(c|e))d" == add_branch_or("a(b|c)d", "e", idx = 5)
        @test "a(b|c)(d|e)" == add_branch_or("a(b|c)d", "e", idx = 7)

        @test "((a|e)b|c)" == add_branch_or("(ab|c)", "e", idx = 2)
        @test "(a(b|e)|c)" == add_branch_or("(ab|c)", "e", idx = 3)
        @test "(ab|(c|e))" == add_branch_or("(ab|c)", "e", idx = 5)
    end

    @testset "Addition when and branch present" begin
        @test "(a|e)[bc]{2}d" == add_branch_or("a[bc]{2}d", "e", idx = 1)
        @test "a([bc]{2}|e)d" == add_branch_or("a[bc]{2}d", "e", idx = 2)
        @test "a[bc]{2}(d|e)" == add_branch_or("a[bc]{2}d", "e", idx = 9)

        @test "a([bcd]{2}|f)e" == add_branch_or("a[bcd]{2}e", "f", idx = 2)
    end
end

@testset "Remove or branch mutation" begin
    @test "abcd" == remove_branch_or("a(bc)d", idx = 2)
    @test "abcd" == remove_branch_or("a(bc)+d", idx = 2)
    @test "abd" == remove_branch_or("a(b|c)d", idx = 2, take_l_branch = true)
    @test "acd" == remove_branch_or("a(b|c)d", idx = 2, take_l_branch = false)

    @test "a(b|c)e" == remove_branch_or("a((b|c)|d)e", idx = 2, take_l_branch = true)
    @test "ade" == remove_branch_or("a((b|c)|d)e", idx = 2, take_l_branch = false)

    @test "a(b|c)g" == remove_branch_or("a((b|c)|[def]{2})g", idx = 2, take_l_branch = true)
    @test "a[def]{2}g" == remove_branch_or("a((b|c)|[def]{2})g", idx = 2, take_l_branch = false)
end

@testset "Add loop mutation" begin
    @test "a+bcd" == add_loop("abcd", idx = 1, idx2 = 1)
    @test "ab+cd" == add_loop("abcd", idx = 2, idx2 = 2)
    @test "abc+d" == add_loop("abcd", idx = 3, idx2 = 3)
    @test "abcd+" == add_loop("abcd", idx = 4, idx2 = 4)

    @test "(ab)+cd" == add_loop("abcd", idx = 1, idx2 = 2)
    @test "(abc)+d" == add_loop("abcd", idx = 1, idx2 = 3)
    @test "(abcd)+" == add_loop("abcd", idx = 1, idx2 = 4)

    @test "a(bc)+d" == add_loop("abcd", idx = 2, idx2 = 3)
    @test "a(bcd)+" == add_loop("abcd", idx = 2, idx2 = 4)
    @test "ab(cd)+" == add_loop("abcd", idx = 3, idx2 = 4)
end

@testset "Remove loop mutation" begin
    @test "abcd" == remove_loop("a+bcd")
    @test "abcd" == remove_loop("ab+cd")
    @test "abcd" == remove_loop("abc+d")
    @test "abcd" == remove_loop("abcd+")

    @test "(a|b)cd" == remove_loop("(a|b)+cd")
    @test "(ab|c)d" == remove_loop("(ab|c)+d")
    @test "(ab|cd)" == remove_loop("(ab|cd)+")

    @test "abcd" == remove_loop("(ab)+cd")
    @test "abcd" == remove_loop("(abc)+d")
    @test "abcd" == remove_loop("(abcd)+")
    @test "abcd" == remove_loop("a(bc)+d")
    @test "abcd" == remove_loop("a(bcd)+")
    @test "abcd" == remove_loop("ab(cd)+")
end

@testset "Add branch and mutation" begin
    @test "a[bc]{2}" == add_branch_and("abc", "c", idx = 2)
    @test "ab[cd]{2}" == add_branch_and("abd", "cd", idx = 3)
    @test "a[bc]{2}d" == add_branch_and("ab+d", "bc", idx = 2)
end

@testset "Pull out mutation" begin
    @test nothing == pull_out("(a[bc]{2}|d)", 1)
    @test nothing == pull_out("(a[bc]{2}|(a|d))", 1)

    @test nothing == pull_out("(a[bc]{2}|a+)", 1)
    @test nothing == pull_out("(a+[bc]{2}|a)", 1)

    @test nothing == pull_out("(a+[bc]{2}a)", 1)
    @test nothing == pull_out("(a[bc]{2}(a|d))", 1)

    @test "ab" == pull_out("(a|ab)", 1)
    @test "ab" == pull_out("(ab|b)", 1)

    @test "a(b|bc)" == pull_out("(ab|abc)", 1)
    @test "(ab|b)c" == pull_out("(abc|bc)", 1)

    @test "a[bc]{2}" == pull_out("(a[bc]{2}|a)", 1)
    @test "[ab]{2}c" == pull_out("([ab]{2}c|c)", 1)

    @test "a+[bc]{2}" == pull_out("(a+[bc]{2}|a+)", 1)
    @test "[ab]{2}c+" == pull_out("([ab]{2}c+|c+)", 1)

    @test "(a|[bc]{2})d" == pull_out("(ad|[bc]{2}d)", 1)
    @test "a([bc]{2}|d)" == pull_out("(a[bc]{2}|ad)", 1)

    @test "[ab]{2}cd" == pull_out("([ab]{2}|[ab]{2}cd)", 1)
    @test "ab[cd]{2}" == pull_out("(ab[cd]{2}|[cd]{2})", 1)

    # Should we support this?
    # @test "[ab]{2}+cd" == pull_out("([ab]{2}+|[ab]{2}+cd)", 1)
    # @test "ab[cd]{2}+" == pull_out("(ab[cd]{2}|[ab]{2}+)", 1)

    @test "(a|b)cd" == pull_out("((a|b)|(a|b)cd)", 1)
    @test "ab(c|d)" == pull_out("(ab(c|d)|(c|d))", 1)

    @test "(a|b)+cd" == pull_out("((a|b)+|(a|b)+cd)", 1)
    @test "ab(c|d)+" == pull_out("(ab(c|d)+|(c|d)+)", 1)

    @test nothing == pull_out("((a|b)|(a|b)+cd)", 1)
    @test nothing == pull_out("(ab(c|d)|(c|d)+)", 1)

    @test "(ab|[cd]{2})([ef]{2}|gh)" == pull_out("((ab|[cd]{2})[ef]{2}|(ab|[cd]{2})gh)", 1)
    @test "((ab+|[cd]{2})|(ab|[cd]{2})gh)[ef]{2}" == pull_out("((ab+|[cd]{2})[ef]{2}|(ab|[cd]{2})gh[ef]{2})", 1)
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
