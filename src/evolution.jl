module evolution

export add_loop, add_branch, add_state, crossover, init_population,
mutate

include("parameters.jl")

using Random: rand
using StatsBase: sample, Weights

@enum Mutation begin
    ADD_EVENT
    REMOVE_EVENT
    ADD_BRANCH
    REMOVE_BRANCH
    ADD_LOOP
    REMOVE_LOOP
end

const ALL_MUTATIONS = Dict(
    ADD_EVENT => ADD_EVENT_RATE
    REMOVE_EVENT => REMOVE_EVENT_RATE
    ADD_BRANCH => ADD_BRANCH_RATE,
    REMOVE_BRANCH => REMOVE_BRANCH_RATE,
    ADD_LOOP => ADD_LOOP_RATE,
    REMOVE_LOOP => REMOVE_LOOP_RATE,
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

function remove_loop(chromo::String)::String
    plus_ids = findall(l -> l == '+', chromo)
    isempty(plus_ids) && return "+-"

    idx = rand(plus_ids)

    return chromo[1: idx - 1] * chromo[idx + 1: end]
end

function add_branch(chromo::String, letters::String)::String
    letter = rand(letters)
    state_ids = findall(l -> isletter(l) && l != letter, collect(chromo))

    isempty(state_ids) && return "(+"
    idx = rand(state_ids)

    return chromo[1: idx - 1] * '(' * chromo[idx] * '|' * letter * ')' *
        chromo[idx + 1: end]
end

function remove_branch(chromo::String)::String
    bracket_ids = findall(l -> l == '(', chromo)
    isempty(bracket_ids) && return "(-"

    idx = rand(bracket_ids)

    l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)

    left_branch = chromo[l_bracket: pipe-1]
    right_branch = chromo[pipe+1: r_bracket]

    if plus
        r_bracket += 1
    end

    return chromo[1: idx - 1] * rand([left_branch, rigth_branch]) *
        chromo[c_bracket_idx + 1: end]
end

function add_state(chromo::String, letters::String)::String
    idx = rand(1:length(chromo))
    letter = rand(letters)

    return chromo[1: idx] * letter * chromo[idx + 1: end]
end

function remove_event(chromo::String)::String
    event = rand(findall(isletter, chromo))
    r_event = event

    if length(chromo) > idx && chromo[event+1] == '+'
        r_event += 1
    end

    if (event > 1 && chromo[event-1] == '(') ||
        (length(chromo) > r_event && chromo[r_event+1] == '|')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, event)

        if plus
            r_bracket += 1
        end
        return chromo[1:l_bracket-1] * chromo[pipe+1:r_bracket-1] * chromo[r_bracket+1:end]
    elseif (event > 1 && chromo[event-1] == '|') ||
        (length(chromo) > r_event && chromo[r_event+1] == ')')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, event)

        if plus
            r_bracket += 1
        end
        return chromo[1:l_bracket-1] * chromo[l_bracket+1:pipe-1] * chromo[r_bracket+1:end]
    else
        return chromo[1:event-1] * chromo[r_event+1]
    end
end

function crossover(chromo1::String, chromo2::String)::String
    event1 = rand(findall(l -> isletter(l) || l == '(', collect(chromo1)))
    event2 = rand(findall(l -> isletter(l) || l == '(', collect(chromo2)))

    sub_chromo1 = if chromo1[event1] == '('
        _, _, r_bracket, plus = find_brackets(chromo1, event1)
        if plus
            r_bracket += 1
        end
        chromo1[1: event1 - 1] * chromo1[r_bracket + 1: end]
    else
        chromo1[1: event1 - 1] * chromo1[event1 + 1: end]
    end
    event1 -= 1

    sub_chromo2 = if chromo2[event2] == '('
        _, _, r_bracket, plus = find_brackets(chromo2, event2)
        if plus
            r_bracket += 1
        end
        chromo2[event2: r_bracket]
    else
        chromo2[event2]
    end

    return if isempty(sub_chromo1) && isempty(sub_chromo2)
         string(chromo1)
    elseif isempty(sub_chromo1)
         string(sub_chromo2)
    elseif isempty(sub_chromo2)
         string(sub_chromo1)
    else
         sub_chromo1[1: event1] * sub_chromo2 * sub_chromo1[event1 + 1: end]
    end
end

function find_brackets(chromo::String, opening_bracket_idx::Int)
   left_bracket, pipe, right_bracket = nothing, nothing, nothing

    bracket_counter = 0
    for i = idx:length(chromo)
        if chromo[i] == ')' && bracket_counter == 0
            right_bracket = i
            break
        elseif chromo[i] == '|' && bracket_counter == 0
            pipe = i
        elseif chromo[i] == ')'
            bracket_counter -= 1
        elseif chromo[i] == '(' && i != idx
            bracket_counter += 1
        elseif chromo[i] == '+' && i == idx && chromo[i-1] == ')'
            right_bracket = i - 1
            bracket_counter -= 1
            break
        end
    end

    for i = idx:-1:1
        if chromo[i] == '(' && bracket_counter == 0
            left_bracket = i
            break
        elseif chromo[i] == '|' && bracket_counter == 0
            pipe = i
        elseif chromo[i] == '('
            bracket_counter -= 1
        elseif chromo[i] == ')' && i != idx
            bracket_counter += 1
        end
    end

    plus_presence = if right_bracket != nothing
        length(chromo) > right_bracket && chromo[right_bracket+1] == '+'
    else
        nothing
    end

    return left_bracket, pipe, right_bracket, plus_presence
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
            elseif mutation == ADD_EVENT
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
