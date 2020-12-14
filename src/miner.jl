module miner

include("scoring.jl")
include("evolution.jl")
include("parameters.jl")

using Printf: @printf
using StatsBase: sample
using DataStructures: OrderedDict

using .scoring
using .evolution

# parse log file
logs = readlines(LOG_FILE)

# extract all possible letters representing events
letters = join(unique(join(logs)))

counterexamples = init_counterexamples(letters, logs)

old_population = init_population(letters, logs)

top_rank_list = OrderedDict(
    score_population(old_population, logs, counterexamples)
)
sort!(top_rank_list, byvalue = true)

for i = 1:ITERATION_NUMBER
    println("[GENERATION $(i)]")
    global top_rank_list, old_population

    # reproduce best chromosomes
    new_population = collect(keys(top_rank_list))

    # replenish population with mutated chromosomes
    new_population = mutate(new_population, old_population, letters)

    # evaluate each chromosome's fitness
    new_scores = score_population(new_population, logs, counterexamples)

    # update top rank list
    for (chromo, score) in new_scores
        if !(chromo in keys(top_rank_list))
            top_rank_list[chromo] = score
        end
    end
    sort!(top_rank_list, byvalue = true)

    # trim top rank list if needed
    while length(top_rank_list) > REPRODUCTION_SIZE
        pop!(top_rank_list)
    end

    # proceed to next generation
    old_population = new_population

    # print top chromosome
    chromo, penalty = first(top_rank_list)
    @printf("'%s' => (%.3f, %.3f)\n", chromo, penalty...)
end

end
