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
        possible_ids = Vector{Int}()
        for (i, e) in enumerate(chromo)
            if isletter(e) || occursin(e, "[(|)")
                push!(possible_ids, i - 1, i)
            elseif e == '}'
                push!(possible_ids, i)
            elseif e == ']'
                push!(possible_ids, i - 1)
            end
        end
        idx = rand(possible_ids)
    end

    event = rand(events)

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
    new_event = rand(events)

    possible_ids = Vector{Int}()
    inside_and = false
    for (i, e) in enumerate(chromo)
        if !inside_and && isletter(e) && e != new_event
            push!(possible_ids, i)
        elseif e == '['
            inside_and = true
        elseif e == '}'
            inside_and = false
        end
    end

    isempty(possible_ids) && return nothing

    if idx == nothing
        idx = rand(possible_ids)
    end

    return chromo[1:idx-1] * '(' * chromo[idx] * '|' * new_event * ')' * chromo[idx+1:end]
end

function remove_branch(chromo::String; idx = nothing)::Union{String,Nothing}
    bracket_ids = findall(l -> l == '(', chromo)
    isempty(bracket_ids) && return nothing

    if idx == nothing
        idx = rand(bracket_ids)
    end

    l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)

    left_branch = chromo[l_bracket+1:pipe-1]
    right_branch = chromo[pipe+1:r_bracket-1]

    r_shift = plus ? 1 : 0
    return chromo[1:l_bracket-1] *
           rand([left_branch, right_branch]) *
           chromo[r_bracket+r_shift:end]
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

function add_and(chromo::String, events::String; idx = nothing)::Union{String,Nothing}

    if idx == nothing
        possible_ids = Vector{Int}()
        inside_and = false
        for (i, e) in enumerate(chromo)
            if isletter(e) && !inside_and
                push!(possible_ids, i)
            elseif e == '['
                inside_and = true
            elseif e == '}'
                inside_and = false
            end
        end
        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end

    new_sub_chromo, r_idx =
        if idx < length(chromo) && isletter(chromo[idx+1]) && chromo[idx+1] != chromo[idx]
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

function crossover(
    chromo1::String,
    chromo2::String;
    idx1 = nothing,
    idx2 = nothing,
)::Union{String,Nothing}
    if idx1 == nothing
        possible_ids = Vector{Int}()
        inside_and = false
        for (i, e) in enumerate(chromo1)
            if !inside_and && (isletter(e) || e == '(')
                push!(possible_ids, i)
            elseif e == '['
                push!(possible_ids, i)
                inside_and = true
            elseif e == '}'
                inside_and = false
            end
        end
        isempty(possible_ids) && return nothing
        idx1 = rand(possible_ids)
    end
    if idx2 == nothing
        possible_ids = Vector{Int}()
        inside_and = false
        for (i, e) in enumerate(chromo2)
            if !inside_and && (isletter(e) || e == '(')
                push!(possible_ids, i)
            elseif e == '['
                push!(possible_ids, i)
                inside_and = true
            elseif e == '}'
                inside_and = false
            end
        end
        isempty(possible_ids) && return nothing
        idx2 = rand(possible_ids)
    end

    sub_chromo1 = if chromo1[idx1] == '('
        _, _, r_bracket, plus = find_brackets(chromo1, idx1)
        if plus
            r_bracket += 1
        end
        chromo1[r_bracket+1:end]
    elseif chromo1[idx1] == '['
        r_and = findfirst(e -> e == '}', chromo1[idx1+6:end]) + length(chromo1[1:idx1+5])
        chromo1[r_and+1:end]
    else
        r_shift = idx1 < length(chromo1) && chromo1[idx1+1] == '+' ? 2 : 1
        chromo1[idx1+r_shift:end]
    end
    sub_chromo1 = chromo1[1:idx1-1] * sub_chromo1
    idx1 -= 1

    sub_chromo2 = if chromo2[idx2] == '('
        _, _, r_bracket, plus = find_brackets(chromo2, idx2)
        if plus
            r_bracket += 1
        end
        chromo2[idx2:r_bracket]
    elseif chromo2[idx2] == '['
        r_and = findfirst(e -> e == '}', chromo2[idx2+6:end]) + length(chromo2[1:idx2+5])
        chromo2[idx2:r_and]
    else
        r_shift = idx2 < length(chromo2) && chromo2[idx2+1] == '+' ? 1 : 0
        chromo2[idx2:idx2+r_shift]
    end

    return if isempty(sub_chromo1)
        string(sub_chromo2)
    else
        sub_chromo1[1:idx1] * sub_chromo2 * sub_chromo1[idx1+1:end]
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

    plus_present = if right_bracket != nothing
        right_bracket < length(chromo) && chromo[right_bracket+1] == '+'
    else
        nothing
    end

    return left_bracket, pipe, right_bracket, plus_present
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
                @debug "Chromo too small: '$(mutation == nothing ? "crossover" : mutation)': '$(chromo)' -> '$(new_chromo)'"
            end
        catch e
            @debug "Faulty chromo: '$(mutation)': '$(chromo)' -> '$(new_chromo)'"
        end
    end

    return new_population
end

end
