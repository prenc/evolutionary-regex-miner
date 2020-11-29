module miner

include("scoring.jl")
include("evolution.jl")
include("parameters.jl")

using StatsBase: sample
using DataStructures: OrderedDict
using .scoring
using .evolution

lines = readlines(LOG_FILE)

println("Initial regex: ", EXAMPLE_REGEX)
println("Score: ", score(EXAMPLE_REGEX, lines))

println()
println("add loop:")
for _ in 1:5
    println(add_loop(EXAMPLE_REGEX))
end

println()
println("add branch:")

for _ in 1:5
    println(add_branch(EXAMPLE_REGEX, LETTERS))
end

println()
println("add state:")

for _ in 1:5
    println(add_state(EXAMPLE_REGEX, LETTERS))
end

println()
println("crossover:")

for _ in 1:5
    println(crossover(EXAMPLE_REGEX, "ab(d(a|g)|(a|g)c)+"))
end

println()
# algorithm
println("Actual algorithm")

old_population = init_population(LETTERS, POPULATION_SIZE)

# todo change to SortedDict which sorts by value
top_rank_list = OrderedDict{String,Float64}(score_population(old_population, lines))

for i in 1:ITERATION_NUMBER
    println("[ITERATION $(i)]")
    global top_rank_list, old_population

    new_population = Vector{String}()
    println("Top 5 chromosomes:")
    for (i, (chromo, penalty)) in enumerate(pairs(top_rank_list))
        if i <= 5
            println("chromo: ", chromo, ", penalty: ", penalty)
        end
        if length(new_population) <= REPRODUCTION_SIZE
            push!(new_population, chromo)
        else
            break
        end
    end

    new_population = mutate(new_population, old_population)

    new_scores = score_population(new_population, lines)

    top_rank_list = OrderedDict(
        sort(
            vcat(collect(pairs(top_rank_list)), new_scores),
            by = x -> x[2]
       )
    )

    while length(top_rank_list) >= TOP_LIST_SIZE
        pop!(top_rank_list)
    end

    println("New chromosomes sample:")
    for (chromo, penalty) in sample(new_scores, 5)
        println("chromo: ", chromo, ", penalty: ", penalty)
    end

    old_population = new_population
end
end
