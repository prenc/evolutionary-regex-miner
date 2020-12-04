module miner

include("scoring.jl")
include("evolution.jl")
include("parameters.jl")

using StatsBase: sample
using DataStructures: OrderedDict

using .scoring
using .evolution

logs = readlines(LOG_FILE)
letters = join(unique(join(logs)))

old_population = init_population(letters, POPULATION_SIZE)

top_rank_list = OrderedDict{String,Float64}(score_population(old_population, logs))
sort!(top_rank_list, byvalue = true)

for i in 1:ITERATION_NUMBER
    println("[GENERATION $(i)]")
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

    new_population = mutate(new_population, old_population, letters)

    new_scores = score_population(new_population, logs)

    for (chromo, score) in new_scores
        if !(chromo in keys(top_rank_list))
            top_rank_list[chromo] = score
        end
    end
    sort!(top_rank_list, byvalue = true)

    while length(top_rank_list) > TOP_LIST_SIZE
        pop!(top_rank_list)
    end

    println("New chromosomes sample:")
    for (chromo, penalty) in sample(new_scores, 5)
        println("chromo: ", chromo, ", penalty: ", penalty)
    end

    old_population = new_population
end

end
