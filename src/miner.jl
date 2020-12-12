module miner

include("scoring.jl")
include("evolution.jl")
include("parameters.jl")

using Printf: @printf
using StatsBase: sample
using DataStructures: OrderedDict

using .scoring
using .evolution

logs = readlines(LOG_FILE)
letters = join(unique(join(logs)))

counterexamples = init_counterexamples(letters, logs)

old_population = init_population(letters, POPULATION_SIZE, logs)

top_rank_list = OrderedDict(score_population(old_population, logs, counterexamples))
sort!(top_rank_list, byvalue = true)

for i in 1:ITERATION_NUMBER
    println("[GENERATION $(i)]")
    global top_rank_list, old_population

    new_population = Vector{String}()

    for (i, (chromo, penalty)) in enumerate(pairs(top_rank_list))
        if length(new_population) <= REPRODUCTION_SIZE
            push!(new_population, chromo)
        else
            break
        end
    end

    new_population = mutate(new_population, old_population, letters)

    new_scores = score_population(new_population, logs, counterexamples)

    for (chromo, score) in new_scores
        if !(chromo in keys(top_rank_list))
            top_rank_list[chromo] = score
        end
    end
    sort!(top_rank_list, byvalue = true)

    while length(top_rank_list) > TOP_LIST_SIZE
        pop!(top_rank_list)
    end

    old_population = new_population

    for (i, (chromo, penalty)) in enumerate(pairs(top_rank_list))
        if i <= 1
            @printf("'%s' => (%.3f, %.3f)\n", chromo, penalty...)
        else
            break
        end
    end
end

end
