module evolution

export init_population, mutate

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
    ADD_AND
    REMOVE_AND
end

const ALL_MUTATIONS = Dict(
    ADD_EVENT => ADD_EVENT_RATE,
    REMOVE_EVENT => REMOVE_EVENT_RATE,
    ADD_BRANCH => ADD_BRANCH_RATE,
    REMOVE_BRANCH => REMOVE_BRANCH_RATE,
    ADD_LOOP => ADD_LOOP_RATE,
    REMOVE_LOOP => REMOVE_LOOP_RATE,
    ADD_AND => ADD_AND_RATE,
    REMOVE_AND => REMOVE_AND_RATE,
)

function add_event(chromo::String, events::String; idx = nothing)::String
    if idx == nothing
        possible_ids = findall(e -> !isdigit(e) && !occursin(e, "]{"), chromo)
        idx = rand(possible_ids)
    end
    event = rand(events)

    idx -= rand([0, 1])

    return chromo[1:idx] * event * chromo[idx+1:end]
end

function remove_event(chromo::String; idx = nothing)::String
    event = idx == nothing ? rand(findall(isletter, chromo)) : idx
    r_event = event

    if event < length(chromo) && chromo[event+1] == '+'
        r_event += 1
    end

    return if (event > 1 && chromo[event-1] == '(') ||
              (r_event < length(chromo) && chromo[r_event+1] == '|')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, event)
        r_shift = plus ? 2 : 1

        chromo[1:l_bracket-1] * chromo[pipe+1:r_bracket-1] * chromo[r_bracket+r_shift:end]
    elseif (event > 1 && chromo[event-1] == '|') ||
           (r_event < length(chromo) && chromo[r_event+1] == ')')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, event)
        r_shift = plus ? 2 : 1

        chromo[1:l_bracket-1] * chromo[l_bracket+1:pipe-1] * chromo[r_bracket+r_shift:end]
    elseif event > 1 && chromo[event-1] == '['
        r_bracket = findfirst(e -> e == ']', chromo[event+1:end]) + length(chromo[1:event])
        new_chromo = chromo[1:event-1] * chromo[event+1:end]

        if r_bracket - event <= 2
            remove_and(new_chromo, idx = event - 1)
        else
            new_chromo
        end
    elseif event < length(chromo) && chromo[event+1] == ']'
        l_bracket =
            length(chromo[1:event-1]) - findfirst(e -> e == '[', reverse(chromo[1:event-2]))
        new_chromo = chromo[1:event-1] * chromo[event+1:end]

        if event - l_bracket <= 2
            remove_and(new_chromo, idx = l_bracket)
        else
            new_chromo
        end
    else
        chromo[1:event-1] * chromo[r_event+1:end]
    end
end

function add_branch(chromo::String, events::String; idx = nothing)::Union{String,Nothing}
    event = rand(events)

    event_ids = Vector{Int}()
    can_add = true
    for i in 1:length(chromo)
        if can_add && isletter(chromo[i]) && chromo[i] != event
            push!(event_ids, i)
        elseif chromo[i] == '['
            can_add = false
        elseif chromo[i] == '}'
            can_add = true
        end
    end

    isempty(event_ids) && return nothing
    if idx == nothing
        idx = rand(event_ids)
    end

    return chromo[1:idx-1] * '(' * chromo[idx] * '|' * event * ')' * chromo[idx+1:end]
end

function remove_branch(chromo::String; idx = nothing)::Union{String,Nothing}
    bracket_ids = findall(l -> l == '(', chromo)
    isempty(bracket_ids) && return nothing

    if idx == nothing
        idx = rand(bracket_ids)
    end

    l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)

    left_branch = chromo[l_bracket:pipe-1]
    right_branch = chromo[pipe+1:r_bracket]

    return if plus
        chromo[1:idx] * rand([left_branch, right_branch]) * chromo[r_bracket:end]
    else
        chromo[1:idx-1] * rand([left_branch, right_branch]) * chromo[r_bracket+1:end]
    end
end

function add_loop(chromo::String; idx = nothing)::Union{String,Nothing}
    possible_ids = Vector()
    for i = 1:length(chromo)
        if (isletter(chromo[i]) || chromo[i] == ')') &&
           (length(chromo) == i || chromo[i+1] != '+')
            push!(possible_ids, i)
        end
    end

    isempty(possible_ids) && return nothing
    if idx == nothing
        idx = rand(possible_ids)
    end

    return chromo[1:idx] * '+' * chromo[idx+1:end]
end

function remove_loop(chromo::String; idx = nothing)::Union{String,Nothing}
    plus_ids = findall(l -> l == '+', chromo)
    isempty(plus_ids) && return nothing

    if idx == nothing
        idx = rand(plus_ids)
    end

    return chromo[1:idx-1] * chromo[idx+1:end]
end

function add_and(chromo::String, events::String; idx = nothing)::String
    if idx == nothing
        event_ids = findall(isletter, chromo)
        idx = rand(event_ids)
    end

    new_sub_chromo, r_idx = if idx < length(chromo) && isletter(chromo[idx+1])
        chromo[idx:idx+1], 2
    else
        new_events = Set(collect(events))
        pop!(new_events, chromo[idx])
        chromo[idx] * rand(new_events), 1
    end

    return chromo[1:idx-1] * '[' * new_sub_chromo * "]{2}" * chromo[idx+r_idx:end]
end

function remove_and(chromo::String; idx = nothing)::Union{String,Nothing}
    bracket_ids = findall(e -> e == '[', chromo)
    isempty(bracket_ids) && return nothing

    if idx == nothing
        idx = rand(bracket_ids)
    end
    r_bracket = findfirst(e -> e == ']', chromo[idx+2:end]) + length(chromo[1:idx+1])

    return chromo[1:idx-1] * chromo[idx+1:r_bracket-1] * chromo[r_bracket+4:end]
end

function crossover(chromo1::String, chromo2::String)::Union{String,Nothing}
    event1_idx = rand(findall(e -> isletter(e) || e == '(', collect(chromo1)))
    event2_idx = rand(findall(e -> isletter(e) || e == '(', collect(chromo2)))

    sub_chromo1 = if chromo1[event1_idx] == '('
        _, _, r_bracket, plus = find_brackets(chromo1, event1_idx)
        if plus
            r_bracket += 1
        end
        chromo1[r_bracket+1:end]
    else
        chromo1[event1_idx+1:end]
    end
    sub_chromo1 = chromo1[1:event1_idx-1] * sub_chromo1
    event1_idx -= 1

    sub_chromo2 = if chromo2[event2_idx] == '('
        _, _, r_bracket, plus = find_brackets(chromo2, event2_idx)
        if plus
            r_bracket += 1
        end
        chromo2[event2_idx:r_bracket]
    else
        chromo2[event2_idx]
    end

    return if isempty(sub_chromo1)
        string(sub_chromo2)
    else
        sub_chromo1[1:event1_idx] * sub_chromo2 * sub_chromo1[event1_idx+1:end]
    end
end

function find_brackets(chromo::String, idx::Int)
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
        right_bracket < length(chromo) && chromo[right_bracket+1] == '+'
    else
        nothing
    end

    return left_bracket, pipe, right_bracket, plus_presence
end

function init_population(events::String, population_size::Int)::Vector{String}
    population = Vector{String}()

    while length(population) < population_size
        chromo_size = rand(MIN_CHROMO_SIZE:INITIAL_MAX_CHROMO_SIZE)
        push!(population, join(sample(collect(events), chromo_size)))
    end

    return population
end

function mutate(
    new_population::Vector{String},
    old_population::Vector{String},
    events::String,
)::Vector{String}
    while length(new_population) <= POPULATION_SIZE
        chromo = rand(old_population)

        evolution_chance = rand()
        new_chromo = if evolution_chance <= CROSSOVER_RATE
            crossover(chromo, rand(old_population))
        else
            mutation = sample(
                collect(keys(ALL_MUTATIONS)),
                Weights(collect(values(ALL_MUTATIONS))),
            )
            if mutation == ADD_EVENT
                add_event(chromo, events)
            elseif mutation == REMOVE_EVENT
                remove_event(chromo)
            elseif mutation == ADD_BRANCH
                add_branch(chromo, events)
            elseif mutation == REMOVE_BRANCH
                remove_branch(chromo)
            elseif mutation == ADD_LOOP
                add_loop(chromo)
            elseif mutation == REMOVE_LOOP
                remove_loop(chromo)
            elseif mutation == ADD_AND
                add_and(chromo, events)
            elseif mutation == REMOVE_AND
                remove_and(chromo)
            else
                error("Unsupported mutation encounterd: '$(mutation)'")
            end
        end

        if new_chromo == nothing
            @debug "Unable to mutate: '$(mutation)': '$(chromo)'"
            continue
        end
        try
            Regex(new_chromo)
            if length(new_chromo) >= MIN_CHROMO_SIZE
                push!(new_population, new_chromo)
            else
                @debug "Chromo too small: '$(mutation)': '$(chromo)' -> '$(new_chromo)'"
            end
        catch e
            @debug "Faulty chromo: '$(mutation)': '$(chromo)' -> '$(new_chromo)'"
        end
    end

    return new_population
end

end
