module evolution

export add_loop, add_branch, add_state, crossover, init_population,
mutate

include("parameters.jl")

using Random: rand
using StatsBase: sample, Weights

@enum Mutation ADD_BRANCH ADD_LOOP ADD_STATE
const ALL_MUTATIONS = Dict(
    ADD_BRANCH::Mutation => ADD_BRANCH_RATE,
    ADD_LOOP::Mutation => ADD_LOOP_RATE,
    ADD_STATE::Mutation => ADD_STATE_RATE
)

function add_loop(chromosome::String)::String
    possible_ids = findall(l -> isletter(l) || l == ')', collect(chromosome))

    idx = rand(1:length(chromosome))
    return chromosome[1: idx] * '+' * chromosome[idx + 1: length(chromosome)]
end

function add_branch(chromosome::String, letters::String)::String
    letter = rand(letters)
    state_ids = findall(l -> isletter(l) && l != letter, collect(chromosome))
    isempty(state_ids) && return chromosome
    idx = rand(state_ids)

    return chromosome[1: idx - 1] * '(' * chromosome[idx] * '|' * letter * ')' *
        chromosome[idx + 1: length(chromosome)]
end

function add_state(chromosome::String, letters::String)::String
    idx = rand(1:length(chromosome))
    letter = rand(letters)

    return chromosome[1: idx] * letter * chromosome[idx + 1: length(chromosome)]
end

function crossover(chromosome1::String, chromosome2::String)::String
    state1 = rand(findall(l -> isletter(l) || l == '(', collect(chromosome1)))
    state2 = rand(findall(l -> isletter(l) || l == '(', collect(chromosome2)))

    sub_chromo1 = if chromosome1[state1] == '('
        c_state1 = find_closing_bracket(chromosome1, state1)
        chromosome1[1: state1 - 1] * chromosome1[c_state1 + 1: end]
    else
        new_chrome = chromosome1[1: state1 - 1] * chromosome1[state1 + 1: end]
        new_chrome
    end
    state1 -= 1

    sub_chromo2 = if chromosome2[state2] == '('
        c_state2 = find_closing_bracket(chromosome2, state2)
        chromosome2[state2: c_state2]
    else
        chromosome2[state2]
    end

    return if isempty(sub_chromo1) && isempty(sub_chromo2)
         string(chromosome1)
    elseif isempty(sub_chromo1)
         string(sub_chromo2)
    elseif isempty(sub_chromo2)
         string(sub_chromo1)
    else
         sub_chromo1[1: state1] * sub_chromo2 * sub_chromo1[state1 + 1: end]
    end
end

function find_closing_bracket(chromosome::String, opening_bracket_idx::Int)::Union{Int,Nothing}
    closing_bracket_idx = nothing
    opening_bracket_counter = 0
    for i in opening_bracket_idx + 1: length(chromosome)
        if chromosome[i] == ')' && opening_bracket_counter == 0
            closing_bracket_idx = i
            break
        elseif chromosome[i] == ')'
            opening_bracket_counter -= 1
        elseif chromosome[i] == '('
            opening_bracket_counter += 1
        end
    end
    if closing_bracket_idx < length(chromosome) && chromosome[closing_bracket_idx + 1] == '+'
        closing_bracket_idx += 1
    end
    return closing_bracket_idx
end

function init_population(letters::String, population_size::Int)::Vector{String}
    population = Vector{String}()

    while length(population) < population_size
        chromo_size = rand(MIN_CHROMO_SIZE:MAX_CHROMO_SIZE)
        push!(population, join(sample(collect(letters), chromo_size)))
    end

    return population
end

function mutate(new_population::Vector{String}, old_population::Vector{String})::Vector{String}
    while length(new_population) <= POPULATION_SIZE
        chromo = rand(old_population)
        evolution_chance = rand()
        new_chromo = if evolution_chance <= CROSSOVER_RATE
            crossover(chromo, rand(old_population))
        else
            mutation = sample(
                collect(keys(ALL_MUTATIONS)),
                Weights(collect(values(ALL_MUTATIONS)))
            )
            if mutation == ADD_LOOP
                add_loop(chromo)
            elseif mutation == ADD_BRANCH
                add_branch(chromo, LETTERS)
            elseif mutation == ADD_STATE
                add_state(chromo, LETTERS)
            else
                error("Unsupported mutation encounterd! $(mutation)")
            end
        end

        try
            Regex(new_chromo)
            if length(new_chromo) >= MIN_CHROMO_SIZE
                push!(new_population, new_chromo)
            end
        catch e
        end
    end
    return new_population
end

end
