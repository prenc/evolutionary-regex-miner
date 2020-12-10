module scoring

export score, score_population

include("parameters.jl")

using Statistics
using Base.Threads: @spawn, fetch

function score(regex::String, logs::Vector{String})
    log_fitness = Vector()
    fitness = 1
    for log in logs
        try
            matches = collect(eachmatch(Regex(regex), log))
            push!(
                log_fitness,
                isempty(matches) ? 1 :
                1 - findmax(map(m -> length(m.match), matches))[1] / length(log),
            )
            fitness = mean(log_fitness)
        catch e
            @debug "Chromosome too complex: '$(regex)'"
        end
    end

    event_number = 0
    event_inside_and = 0

    branch_or_number = sum(get_bracket_levels(regex))
    branch_and_number = 0
    loop_number = 0

    inside_and = false
    for c in regex
        if c == '+'
            loop_number += 1
        elseif isletter(c)
            if inside_and
                event_inside_and += 1
            else
                event_number += 1
            end
        elseif c == '['
            branch_and_number += 1
            inside_and = true
        elseif c == '}'
            inside_and = false
        end
    end

    precision =
        event_inside_and * EVENT_INSIDE_AND_PENALTY +
        branch_or_number * BRANCH_PENALTY +
        branch_and_number * AND_PENALTY +
        loop_number * LOOP_PENALTY

    simplicity = event_number * EVENT_PENALTY

    @debug "Score '$(regex)': fitness: '$(fitness)'," *
           " precision: '$(precision)', simplicity: '$(simplicity)'"
    return regex, (fitness, precision + simplicity)
end

function score_population(population::Vector{String}, logs::Vector{String})
    scored_population = map(fetch, map(chromo -> @spawn(score(chromo, logs)), population))

    precisions = [score for (_, (_, score)) in scored_population]
    max_precision = max(precisions...)

    return [(chromo, fitness + precision / max_precision / 2) for (chromo, (fitness, precision)) in scored_population]
end

function get_bracket_levels(chromo::String)::Vector{Int}
    bracket_levels = Vector{Int}()
    level_counter = 1
    for e in chromo
        if e == '('
            push!(bracket_levels, level_counter)
            level_counter += 1
        elseif e == ')'
            level_counter -= 1
        end
    end

    return bracket_levels
end

end
