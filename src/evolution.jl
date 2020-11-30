module evolution

export add_loop, add_branch, add_state, crossover, init_population,
mutate

include("parameters.jl")

using Random: rand
using StatsBase: sample, Weights

@enum Mutation ADD_BRANCH ADD_LOOP ADD_STATE REMOVE_BRANCH REMOVE_LOOP
const ALL_MUTATIONS = Dict(
    ADD_BRANCH::Mutation => ADD_BRANCH_RATE,
    REMOVE_BRANCH::Mutation => REMOVE_BRANCH_RATE,
    ADD_LOOP::Mutation => ADD_LOOP_RATE,
    REMOVE_LOOP::Mutation => REMOVE_LOOP_RATE,
    ADD_STATE::Mutation => ADD_STATE_RATE
)

function add_loop(chromo::String)::String
    possible_ids = Vector()
    for i in 1:length(chromo)
        if (isletter(chromo[i]) || chromo[i] == ')') &&
               (length(chromo) == i || chromo[i + 1] != '+')
            push!(possible_ids, i)
        end
    end

    isempty(possible_ids) && return "++"
    idx = rand(possible_ids)

    return chromo[1: idx] * '+' * chromo[idx + 1: end]
end

function remove_loop(chromosome::String)::String
    plus_ids = findall(l -> l == '+', chromosome)
    isempty(plus_ids) && return "+-"

    idx = rand(plus_ids)

    return chromosome[1: idx - 1] * chromosome[idx + 1: end]
end

function add_branch(chromosome::String, letters::String)::String
    letter = rand(letters)
    state_ids = findall(l -> isletter(l) && l != letter, collect(chromosome))

    isempty(state_ids) && return "(+"
    idx = rand(state_ids)

    return chromosome[1: idx - 1] * '(' * chromosome[idx] * '|' * letter * ')' *
        chromosome[idx + 1: end]
end

function remove_branch(chromosome::String)::String
    bracket_ids = findall(l -> l == '(', chromosome)
    isempty(bracket_ids) && return "(-"

    idx = rand(bracket_ids)

    c_bracket_idx = find_closing_bracket(chromosome, idx)
    pipe_idx = find_closing_bracket(chromosome, idx, '|')

    left_branch = chromosome[idx + 1: pipe_idx - 1]

    rigth_branch = if chromosome[c_bracket_idx] == '+'
        chromosome[pipe_idx + 1: c_bracket_idx - 2]
    else
        chromosome[pipe_idx + 1: c_bracket_idx - 1]
    end

    return chromosome[1: idx - 1] * rand([left_branch, rigth_branch]) *
        chromosome[c_bracket_idx + 1: end]
end

function add_state(chromosome::String, letters::String)::String
    idx = rand(1:length(chromosome))
    letter = rand(letters)

    return chromosome[1: idx] * letter * chromosome[idx + 1: end]
end

function crossover(chromosome1::String, chromosome2::String)::String
    state1 = rand(findall(l -> isletter(l) || l == '(', collect(chromosome1)))
    state2 = rand(findall(l -> isletter(l) || l == '(', collect(chromosome2)))

    sub_chromo1 = if chromosome1[state1] == '('
        c_state1 = find_closing_bracket(chromosome1, state1)
        chromosome1[1: state1 - 1] * chromosome1[c_state1 + 1: end]
    else
        chromosome1[1: state1 - 1] * chromosome1[state1 + 1: end]
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

function find_closing_bracket(chromosome::String, opening_bracket_idx::Int, c::Char=')')::Union{Int,Nothing}
    closing_bracket_idx = nothing
    opening_bracket_counter = 0
    for i in opening_bracket_idx + 1: length(chromosome)
        if chromosome[i] == c && opening_bracket_counter == 0
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
            elseif mutation == REMOVE_LOOP
                remove_loop(chromo)
            elseif mutation == ADD_BRANCH
                add_branch(chromo, LETTERS)
            elseif mutation == REMOVE_BRANCH
                remove_branch(chromo)
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
            else
                @debug "Too small chromo: '$(mutation)': '$(chromo)' -> '$(new_chromo)'"
            end
        catch e
            @debug "Faulty chromo: '$(mutation)': '$(chromo)' -> '$(new_chromo)'"
        end
    end
    return new_population
end

end
