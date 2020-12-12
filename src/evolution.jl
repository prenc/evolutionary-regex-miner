module evolution

export init_population, mutate

include("parameters.jl")

using Random: rand, shuffle
using StatsBase: sample, Weights

@enum Mutation begin
    ADD_EVENT
    REMOVE_EVENT
    ADD_BRANCH_OR
    REMOVE_BRANCH_OR
    ADD_BRANCH_AND
    REMOVE_BRANCH_AND
    ADD_LOOP
    REMOVE_LOOP
    PULL_OUT
end

const ALL_MUTATIONS = Dict(
    ADD_EVENT => ADD_EVENT_RATE,
    REMOVE_EVENT => REMOVE_EVENT_RATE,
    ADD_BRANCH_OR => ADD_BRANCH_OR_RATE,
    REMOVE_BRANCH_OR => REMOVE_BRANCH_OR_RATE,
    ADD_BRANCH_AND => ADD_BRANCH_AND_RATE,
    REMOVE_BRANCH_AND => REMOVE_BRANCH_AND_RATE,
    ADD_LOOP => ADD_LOOP_RATE,
    REMOVE_LOOP => REMOVE_LOOP_RATE,
    PULL_OUT => PULL_OUT_RATE,
)

function add_event(chromo::String, events::String; idx = nothing)::Union{String,Nothing}
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

    bracket_idx = findnext(e -> occursin(e, "[]"), chromo, idx + 1)
    return if bracket_idx != nothing && chromo[bracket_idx] == ']'
        l_bracket = findprev(e -> e == '[', chromo, idx)
        and_content = chromo[l_bracket+1:bracket_idx-1]
        occursin(event, and_content) && return nothing

        chromo[1:l_bracket] *
        join(sort(collect(and_content * event))) *
        chromo[bracket_idx:end]
    else
        chromo[1:idx] * event * chromo[idx+1:end]
    end
end

function remove_event(chromo::String; idx = nothing)::String
    idx = idx == nothing ? rand(findall(isletter, chromo)) : idx
    r_idx = idx

    if idx < length(chromo) && chromo[idx+1] == '+'
        r_idx += 1
    end

    return if (idx > 1 && chromo[idx-1] == '(') ||
              (r_idx < length(chromo) && chromo[r_idx+1] == '|')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)
        r_shift = plus ? 2 : 1
        if pipe == nothing
            if r_bracket - l_bracket == 3
                chromo[1:l_bracket-1] * chromo[idx+1] * chromo[idx+3:end]
            else
                chromo[1:l_bracket] * chromo[idx+1:end]
            end
        else
            chromo[1:l_bracket-1] * chromo[pipe+1:r_bracket-1] * chromo[r_bracket+r_shift:end]
        end
    elseif (idx > 1 && chromo[idx-1] == '|') ||
           (r_idx < length(chromo) && chromo[r_idx+1] == ')')
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)
        r_shift = plus ? 2 : 1
        if pipe == nothing
            if r_bracket - l_bracket == 3
                chromo[1:l_bracket-1] * chromo[idx-1] * chromo[idx+2:end]
            else
                chromo[1:l_bracket] * chromo[idx+1:end]
            end
        else
            chromo[1:l_bracket-1] * chromo[l_bracket+1:pipe-1] * chromo[r_bracket+r_shift:end]
        end
    elseif idx > 1 && chromo[idx-1] == '['
        r_bracket = findnext(e -> e == ']', chromo, idx + 1)
        new_chromo = chromo[1:idx-1] * chromo[idx+1:end]

        if r_bracket - idx <= 2
            remove_branch_and(new_chromo, idx = idx - 1)
        else
            new_chromo
        end
    elseif idx < length(chromo) && chromo[idx+1] == ']'
        l_bracket = findprev(e -> e == '[', chromo, idx - 2)
        new_chromo = chromo[1:idx-1] * chromo[idx+1:end]

        if idx - l_bracket <= 2
            remove_branch_and(new_chromo, idx = l_bracket)
        else
            new_chromo
        end
    else
        chromo[1:idx-1] * chromo[r_idx+1:end]
    end
end

function add_branch_or(chromo::String, events::String; idx = nothing)::Union{String,Nothing}
    new_event = rand(events)

    if idx == nothing
        possible_ids = Vector{Int}()
        inside_and = false
        for (i, e) in enumerate(chromo)
            if !inside_and && isletter(e) && e != new_event
                push!(possible_ids, i)
            elseif e == '['
                push!(possible_ids, i)
                inside_and = true
            elseif e == '}'
                inside_and = false
            end
        end

        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end
    sub_chromo =  if chromo[idx] == '['
        and_end = findnext(e -> e == '}', chromo, idx+6)
        chromo[idx:and_end] * '|' * new_event * ')' * chromo[and_end+1:end]
    else
        chromo[idx] * '|' * new_event * ')' * chromo[idx+1:end]
    end
        return chromo[1:idx-1] * '(' * sub_chromo
end

function remove_branch_or(chromo::String; idx = nothing, take_l_branch = nothing)::Union{String,Nothing}
    if idx == nothing
        possible_ids = findall(l -> l == '(', chromo)

        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end

    l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)

    r_shift = plus ? 2 : 1
    branch = if pipe != nothing
        left_branch = chromo[l_bracket+1:pipe-1]
        right_branch = chromo[pipe+1:r_bracket-1]
        if left_branch == nothing
            rand([left_branch, right_branch])
        elseif take_l_branch
            left_branch
        else
            right_branch
        end
    else
        chromo[l_bracket+1:r_bracket-1]
    end

    return chromo[1:l_bracket-1] * branch * chromo[r_bracket+r_shift:end]
end

function add_loop(chromo::String; idx = nothing, idx2 = nothing)::Union{String,Nothing}
    if idx == nothing
        possible_ids = findall(isletter, chromo)

        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end

    last_letter_idx = idx
    for i in idx+1:length(chromo)
        if isletter(chromo[i])
            last_letter_idx += 1
        else
            break
        end
    end

    random_idx = idx2 == nothing ? rand(idx:last_letter_idx) : idx2
    return if last_letter_idx == idx || random_idx == idx
        chromo[1:idx] * '+' * chromo[idx+1:end]
    else
        chromo[1:idx-1] * '(' * chromo[idx:random_idx] * ")+" * chromo[random_idx+1:end]
    end
end

function remove_loop(chromo::String; idx = nothing)::Union{String,Nothing}
    plus_ids = findall(l -> l == '+', chromo)

    isempty(plus_ids) && return nothing

    if idx == nothing
        idx = rand(plus_ids)
    end

    if chromo[idx-1] == ')'
        l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)
        if pipe == nothing
            return chromo[1:l_bracket-1] * chromo[l_bracket+1:r_bracket-1] * chromo[idx+1:end]
        end
    end
        return chromo[1:idx-1] * chromo[idx+1:end]
end

function add_branch_and(
    chromo::String,
    events::String;
    idx = nothing,
)::Union{String,Nothing}
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
        if idx < length(chromo) && isletter(chromo[idx+1]) && chromo[idx] != chromo[idx+1]
            join(sort(collect(chromo[idx:idx+1]))), 2
        else
            new_events = Set(collect(events))
            pop!(new_events, chromo[idx])
            join(sort([chromo[idx],  rand(new_events)])), 1
        end

        if idx+r_idx <= length(chromo) && chromo[idx+r_idx] == '+'
            r_idx += 1
        end

    return chromo[1:idx-1] * '[' * new_sub_chromo * "]{2}" * chromo[idx+r_idx:end]
end

function remove_branch_and(chromo::String; idx = nothing)::Union{String,Nothing}
    if idx == nothing
        possible_ids = findall(e -> e == '[', chromo)

        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end

    r_bracket = findfirst(e -> e == ']', chromo[idx+2:end]) + length(chromo[1:idx+1])

    return chromo[1:idx-1] *
           join(shuffle(collect(chromo[idx+1:r_bracket-1]))) *
           chromo[r_bracket+4:end]
end

function pull_out(chromo::String, idx = nothing)::Union{String,Nothing}
    if idx == nothing
        possible_ids = findall(e -> e == '(', chromo)

        isempty(possible_ids) && return nothing
        idx = rand(possible_ids)
    end

    l_bracket, pipe, r_bracket, plus = find_brackets(chromo, idx)
    pipe == nothing && return
    r_shift = plus ? 2 : 1

    if chromo[l_bracket+1] == chromo[pipe+1]
        sub_chromo = if isletter(chromo[l_bracket+1])
            plus1 = chromo[l_bracket+2] == '+'

            if (chromo[pipe+2] == '+') == plus1
                r_shift1 = plus1 ? 2 : 1
                if l_bracket + 1 == pipe - r_shift1 || pipe + 1 == r_bracket - r_shift1 # remove or branch
                    chromo[l_bracket+1:l_bracket+r_shift1] * chromo[l_bracket+r_shift1+1:pipe-1] *
                    chromo[pipe+r_shift1+1:r_bracket-1] * chromo[r_bracket+r_shift:end]
                else
                    chromo[l_bracket+1:l_bracket+r_shift1] * '(' * chromo[l_bracket+r_shift1+1:pipe] *
                    chromo[pipe+r_shift1+1:end]
                end
            else
                nothing
            end
        elseif chromo[l_bracket+1] == '['
            and_branch_end1 = findnext(e -> e == '}', chromo, l_bracket+6) # l_bracket+6 first index where } can be found
            and_branch_end2 = findnext(e -> e == '}', chromo, pipe+7) # same as above

            if chromo[l_bracket+1:and_branch_end1] == chromo[pipe+1:and_branch_end2]
                if chromo[pipe-1] == '}' || chromo[r_bracket-1] == '}' # remove or branch
                    chromo[l_bracket+1:and_branch_end1] * chromo[and_branch_end1+1:pipe-1] *
                    chromo[and_branch_end2+1:r_bracket-1] * chromo[r_bracket+r_shift:end]
                else
                    chromo[l_bracket+1:and_branch_end1] * '(' * chromo[and_branch_end1+1:pipe] *
                    chromo[and_branch_end2+1:end]
                end
            else
                nothing
            end
        elseif chromo[l_bracket+1] == '('
            l_bracket1, _, r_bracket1, plus1 = find_brackets(chromo, l_bracket+1)
            l_bracket2, _, r_bracket2, plus2 = find_brackets(chromo, pipe+1)

            if plus1 == plus2 && chromo[l_bracket1:r_bracket1] == chromo[l_bracket2:r_bracket2]
                r_shift1 = plus1 ? 2 : 1
                if chromo[pipe-r_shift1] == ')' || chromo[r_bracket-r_shift1] == ')' # remove or branch
                    chromo[l_bracket1:r_bracket1+r_shift1-1] * chromo[r_bracket1+r_shift1:pipe-1] *
                    chromo[r_bracket2+r_shift1:r_bracket-1] * chromo[r_bracket+r_shift:end]
                else
                    chromo[l_bracket1:r_bracket1] * '(' * chromo[r_bracket1+r_shift1:pipe] *
                    chromo[r_bracket2+r_shift1:end]
                end
            else
                nothing
            end
        else
            nothing
        end
        if sub_chromo != nothing
            return chromo[1:l_bracket-1] * sub_chromo
        end
    end

    if chromo[pipe-1] == chromo[r_bracket-1]
        if chromo[pipe-1] == '+'
            c, plus_sign = chromo[1:pipe-2] * chromo[pipe:r_bracket-2] * chromo[r_bracket:end], '+'
            pipe -= 1
            r_bracket -= 2
        else
            c, plus_sign = chromo, ""
        end
        outer_closing = plus ? ")+" : ")"

        if c[pipe-1] == c[r_bracket-1]
            sub_chromo = if isletter(c[r_bracket-1])
                if r_bracket - pipe == 2 || pipe -  l_bracket == 2 # remove or branch
                    c[1:l_bracket-1] * c[l_bracket+1:pipe-2] *
                    c[pipe+1:r_bracket-1] * plus_sign
                else
                    c[1:pipe-2] * c[pipe:r_bracket-2] * outer_closing *
                    c[r_bracket-1] * plus_sign
                end
            elseif c[r_bracket-1] == '}'
                and_branch_start1 = findprev(e -> e == '[', c, pipe-7) # pipe-6 first index where [ can be found
                and_branch_start2 = findprev(e -> e == '[', c, r_bracket-6) # same as above

                if c[and_branch_start1:pipe-1] == c[and_branch_start2:r_bracket-1]
                    if c[l_bracket+1] == '[' || c[pipe+1] == '[' # remove or branch
                        c[1:l_bracket-1] * c[l_bracket+1:and_branch_start1-1] *
                        c[pipe+1:and_branch_start2-1] * c[and_branch_start2:r_bracket-1]
                    else
                        c[1:and_branch_start1-1] * c[pipe:and_branch_start2-1] * outer_closing *
                        c[and_branch_start2:r_bracket-1]
                    end
                else
                    nothing
                end
            elseif c[r_bracket-1] == ')'
                l_bracket1, _, r_bracket1, plus1 = find_brackets(c, pipe-1)
                l_bracket2, _, r_bracket2, plus2 = find_brackets(c, r_bracket-1)

                if c[l_bracket1:r_bracket1] == c[l_bracket2:r_bracket2]
                    if c[l_bracket+1] == '(' || c[pipe+1] == '(' # remove or branch
                        c[1:l_bracket-1] * c[l_bracket+1:l_bracket1-1] *
                        c[pipe+1:r_bracket2] * plus_sign
                    else
                        c[1:l_bracket1-1] * c[pipe:l_bracket2-1] * outer_closing *
                        c[l_bracket2:r_bracket2] * plus_sign
                    end
                else
                    nothing
                end
            else
                nothing
            end
            if sub_chromo != nothing
                return sub_chromo * c[r_bracket+r_shift:end]
            end
        else
            nothing
        end
    end
    nothing
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

function init_population(
    events::String,
    population_size::Int,
    logs::Vector{String},
)::Vector{String}
    population = Vector{String}()

    while length(population) < population_size
        chromo_size = rand(MIN_CHROMO_SIZE:INITIAL_MAX_CHROMO_SIZE)

        random_log = rand(logs)
        if length(random_log) <= chromo_size
            push!(population, random_log)
        else
            idx = rand(1:length(random_log)-chromo_size)
            push!(population, random_log[idx:idx+chromo_size])
        end
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
            elseif mutation == ADD_BRANCH_OR
                add_branch_or(chromo, events)
            elseif mutation == REMOVE_BRANCH_OR
                remove_branch_or(chromo)
            elseif mutation == ADD_BRANCH_AND
                add_branch_and(chromo, events)
            elseif mutation == REMOVE_BRANCH_AND
                remove_branch_and(chromo)
            elseif mutation == ADD_LOOP
                add_loop(chromo)
            elseif mutation == REMOVE_LOOP
                remove_loop(chromo)
            elseif mutation == PULL_OUT
                pull_out(chromo)
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
