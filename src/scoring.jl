module scoring

include("parameters.jl")

using Statistics

export score, score_population

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
            @debug "Regex failed due to its complexity: '$(regex)'"
        end
    end

    fitness = mean(log_fitness)

    precision = 0 # to be added

    simplicity =
        length(findall(l -> isletter(l) || l == '+', collect(regex))) * STATE_PENALTY

    @debug "Score '$(regex)': fitness: '$(fitness)'," *
           " precision: '$(precision)', simplicity: '$(simplicity)'"
    return fitness + precision + simplicity
end

function score_population(population::Vector{String}, logs::Vector{String})
    return map(chromo -> Pair(chromo, score(chromo, logs)), population)
end

end
