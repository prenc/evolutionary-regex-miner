module miner

include("scoring.jl")
include("evolution.jl")
include("parameters.jl")

using Printf: @printf
using StatsBase: sample
using DataStructures: OrderedDict
using StringDistances: Levenshtein, compare
using Logging

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
    for (new_chromo, new_score) in new_scores
        addition_allowed = true

        for (top_chromo, top_score) in pairs(top_rank_list)
            similarity = compare(new_chromo, top_chromo, Levenshtein())
            if similarity >= SIMILARITY_THRESHOLD
                if new_score < top_score
                    delete!(top_rank_list, top_chromo)
                else
                    addition_allowed = false
                end
                break
            end
        end

        if addition_allowed
            top_rank_list[new_chromo] = new_score
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
    @printf("'%s' => (%.3f)\n", chromo, penalty)

    with_logger(ConsoleLogger(stderr, Logging.Debug)) do
        score(chromo, logs, counterexamples)
    end
end

end
