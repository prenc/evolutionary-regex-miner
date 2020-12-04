module scoring

export score, score_population

include("parameters.jl")

using Statistics

function score(regex::String, logs::Vector{String})::Float64

    log_fitness = Vector()
    for log in logs
        try
            matches = collect(eachmatch(Regex(regex), log))
            push!(
                log_fitness,
                isempty(matches) ? 1 :
                1 - findmax(map(m -> length(m.match), matches))[1] / length(log),
            )
        catch e
            @debug "Chromosome too complex: '$(regex)'"
        end
    end

    fitness = mean(log_fitness)

    event_number = 0
    branch_number = 0
    loop_number = 0

    for idx = 1:length(regex)
        if regex[idx] == '('
            branch_number += 1
        elseif regex[idx] == '+'
            loop_number += 1
        elseif isletter(regex[idx])
            event_number += 1
        end
    end

    precision = branch_number * BRANCH_PENALTY + loop_number * LOOP_PENALTY

    simplicity = event_number * EVENT_PENALTY

    @debug "Score '$(regex)': fitness: '$(fitness)'," *
           " precision: '$(precision)', simplicity: '$(simplicity)'"
    return fitness + precision + simplicity
end

function score_population(population::Vector{String}, logs::Vector{String})
    return map(chromo -> Pair(chromo, score(chromo, logs)), population)
end

end
