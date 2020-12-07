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

old_population = init_population(letters, POPULATION_SIZE, logs)

top_rank_list = OrderedDict{String,Pair{Float64,Float64}}(score_population(old_population, logs))
sort!(top_rank_list, byvalue = true)

for i in 1:ITERATION_NUMBER
    println("[GENERATION $(i)]")
    global top_rank_list, old_population

    new_population = Vector{String}()

    println("Top 5 chromosomes:")
    for (i, (chromo, penalty)) in enumerate(pairs(top_rank_list))
        if i <= 5
            @printf("'%s' => (%.3f, %.5f)\n", chromo, penalty[1], penalty[2])
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

    old_population = new_population
end

end
