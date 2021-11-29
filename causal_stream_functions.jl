using Combinatorics
using PyCall
using Distances
using CSV
using DataFrames

# Structure: output[x][y][z]
# for x, 1 = mechanism, 2 = purview (the mechanism splits have no "reverse twin" splits e.g. [[1],[2]] and [[2],[1]])
# the second level holds all splits of a y sized mechanism/purview
# each z represents one split. This is a 2-element array, where each element is an array specifying the indices for one side of the split of a mechanism/purivew 
struct get_state_map
    n::Int64
    states
    state_index
    split_mechanism
    split_purview

    function get_state_map(n)

        states = Vector[]
        for i in 1:n
            # Make an array that contains arrays of possible system states
            # bitstring converts an int to its binarry version in a long string (lots of zeros)
            # the string is split and we take the last i values and use parse. to convert them to int
            # Reverse the arrays to comply with convention
            push!(states, [reverse(parse.(Int, split(bitstring(x-1), "")[end-(i-1):end])) for x in 1:2^i])
    
        end

         # creating an array with two empty vectors per element
         state_index = [[Vector{Int}(), Vector{Int}()] for i in 1:n]

        # filling up the array with the indices. 
        # The first array for each element represents the system states where that element is 0
        # The second array represent the system states where that element is 1
        for state in 1:2^n, element in 1:n 
            states[n][state][element]==0 ? push!(state_index[element][1], state) : push!(state_index[element][2], state)
        end # nested for loop


        # SPLITS

        split_mechanism = Vector[]
        split_purview = Vector[]

        for i in 1:n

            i == 1 ? push!(split_mechanism, [[[1], Int[]]]) : nothing
            i == 2 ? push!(split_mechanism, [[[1],[2]],[[1,2], Int[]]]) : nothing
    
            if i > 2
    
                split_list = Vector[]
    
                if iseven(i)
                    for x in 1:Int(((i/2)-1))
                        splits = collect(combinations(1:i, x))
                        splits = [[split, setdiff(1:i, split)] for split in splits]
                        for split in splits push!(split_list, split) end
                    end # for x
    
                    half_splits = collect(combinations(1:i, Int(i/2)))
                    n_half = Int(length(half_splits)/2)
    
                    for split in 1:n_half
                        push!(split_list, [half_splits[split], setdiff(1:i, half_splits[split])])
                    end # for split
    
                end # if iseven(i)
    
                
                if isodd(i)
                    for x in 1:Int(((i/2)-.5))
                        splits = collect(combinations(1:i, x))
                        splits = [[split, setdiff(1:i, split)] for split in splits]
                        for split in splits push!(split_list, split) end
                    end # for x
                end # if isodd(i)
                
                push!(split_list, [collect(1:i), Int[]])
    
                push!(split_mechanism, split_list)
    
    
            end # if 1>2
        end # for loop
    
        # Fill out split_purview ... This is where purview splits go
        for i in 1:n
            all_sets = collect(combinations(1:i))
            splits = [[set, setdiff(1:i, set)] for set in all_sets]
            push!(splits, [[],collect(1:i)])
            push!(split_purview, splits)
        end # end for i

        new(n, states, state_index, split_mechanism, split_purview)
    end
end

function convert_TPM(state_map, TPM)

    n = state_map.n
    states = state_map.states[n]
    non_sensor_states = [i[3:n] for i in states]

    TPM = select(TPM, Not(1:3))

    converted_TPM = zeros(2^n,2^n)

    for i in 1:(2^n)
        transfer_state = [TPM[i,col] for col in 1:(n-2)]
        TPM_cols = findall(isequal(transfer_state),non_sensor_states)
        for col in TPM_cols
            converted_TPM[i,col] = 0.25
        end
    end
    
    return converted_TPM
end # function convert_TPM

py"""
from pyemd import emd
import numpy as np

def EMD(d1, d2, metric):
    distance_matrix = np.array(metric)
    results = emd(d1, d2, distance_matrix)
    return(results)
"""

function get_EMD_metric(state_map)

    n = state_map.n

    metric = []
    for i in state_map.states[n]
        push!(metric, convert(Array{Float64,1}, [hamming(collect(i),collect(j)) for j in state_map.states[n]]))
    end

    return metric
end

function get_index_intersection(elements, states, state_index)

    intersection = Vector{Int}()

        @inbounds for i in 1:length(elements)
            element = elements[i]
            state = states[i]
            set = state_index[element][state+1]::Array{Int64,1}
            i==1 ? intersection = set : intersection = intersect(intersection, set)
        end # for i

    return intersection
end

function marginalize_TPM(n, TPM, purview, purview_states, state_index)

    marginalized_TPM = zeros(length(purview_states), 2^n)

    x = 1 
    @inbounds for purview_state in purview_states

        TPM_rows = get_index_intersection(purview, purview_state, state_index)
        
        @inbounds for i in TPM_rows
            marginalized_TPM[x,:] .+= TPM[i,:]
        end # for i

        marginalized_TPM[x,:] ./= length(TPM_rows)

        x += 1
    end # for purview_state

    return marginalized_TPM

end

function cause_rep(TPM, state_map, mechanism, mechanism_state, purview)
    
    n = state_map.n
    system_states = state_map.states[n]::Array{Array{Int64,1},1}
    state_index = state_map.state_index::Array{Array{Array{Int64,1},1},1}


    # get the states to marginalize over
    purview_states = state_map.states[length(purview)]::Array{Array{Int64,1},1}

    marginalized_TPM = marginalize_TPM(n, TPM, purview, purview_states, state_index)

    TPM_columns = get_index_intersection(mechanism, mechanism_state, state_index)

    repetoir = Float64[sum(marginalized_TPM[i,TPM_columns]) for i in 1:length(purview_states)]
    repetoir ./= sum(repetoir) # normalize it into a propability distribution

    return repetoir

end # function cause_rep

function get_split_distance(mech_split_index, pur_split_index, mech1, mech2, mech1_state, mech2_state, pur1, pur2, unconstrained_distributions, TPM, state_map, purview, purview_states, main_rep, EMD_metric, rep_list)
    # SPLIT 1
    if isempty(mech1)
        rep1 = unconstrained_distributions[length(pur1)]
    elseif isempty(pur1)
        nothing # for now
    else
        rep1 = cause_rep(TPM, state_map, mech1 , mech1_state , pur1)
    end

    # SPLIT 2
    if isempty(mech2)
        rep2 = unconstrained_distributions[length(pur2)]
    elseif isempty(pur2)
        nothing # for now
    else
        rep2 = cause_rep(TPM, state_map, mech2 , mech2_state , pur2)
    end
        
    # If one has an empty purview, just use the other.
    if isempty(pur1) 
        combined_rep = rep2
    elseif isempty(pur2) 
        combined_rep = rep1
    else
        # combine destributions
        combined_rep = zeros(2^length(purview))
        j = 1
        for state in purview_states
            state_split1 = state[pur_split_index[1]]
            state_split2 = state[pur_split_index[2]]

            rep1_index = findfirst(isequal(state_split1), state_map.states[length(state_split1)])::Int64
            rep2_index = findfirst(isequal(state_split2), state_map.states[length(state_split2)])::Int64

            combined_rep[j] = rep1[rep1_index] * rep2[rep2_index]
            j += 1

        end # for state
    end # empty / non-empty purview 

    
    if combined_rep in rep_list
        D = -1.0
    else
        D = py"EMD"(main_rep, combined_rep, EMD_metric)::Float64
    end

    return D, combined_rep
end

function get_irreducibility_over_purview(TPM, state_map, mechanism, mechanism_state, purview, EMD_metric, unconstrained_distributions)

    n = state_map.n
    main_rep = cause_rep(TPM, state_map, mechanism, mechanism_state, purview)

    purview_states = state_map.states[length(purview)]::Array{Array{Int64,1},1}

    x = length(mechanism) + length(purview)
    n_distances = 2^(x-1)-1
    distances = zeros(n_distances)

    mech_split_index_list = state_map.split_mechanism[length(mechanism)]
    pur_split_index_list = state_map.split_purview[length(purview)]

    rep_list = []

    i = 1
    @inbounds for mech_split_index::Array{Array{Int64,1},1} in mech_split_index_list, pur_split_index::Array{Array{Int64,1},1} in pur_split_index_list
            
        mech1 = mechanism[mech_split_index[1]]
        mech1_state  = mechanism_state[mech_split_index[1]]
        pur1  = purview[pur_split_index[1]]

        mech2 = mechanism[mech_split_index[2]]
        mech2_state = mechanism_state[mech_split_index[2]]
        pur2  = purview[pur_split_index[2]]

        if !((isempty(mech1)&isempty(pur1)) | (isempty(mech2)&isempty(pur2))) # if no split at all
            
            D, combined_rep = get_split_distance(mech_split_index, pur_split_index, mech1, mech2, mech1_state, mech2_state, pur1, pur2, unconstrained_distributions, TPM, state_map, purview, purview_states, main_rep, EMD_metric, rep_list)
            if D == 0.0 break end
            if D ==-1.0 
                D_index = findfirst(isequal(combined_rep), rep_list)
                distances[i] = distances[D_index]
            else
                distances[i] = D
            end
            i += 1

            push!(rep_list, combined_rep)
        
        end # if no split
    end # for mech_split and pur_split

    irreducibility = min(distances...)

    return irreducibility

end # function get_irreducibility

function get_cause(TPM, state_map, mechanism, mechanism_state, purview_list, EMD_metric, unconstrained_distributions)

    purview_list_length = length(purview_list)
    irreducibility_list = zeros(purview_list_length)

    @inbounds for i in 1:purview_list_length
        purview = purview_list[i]
        irreducibility = get_irreducibility_over_purview(TPM, state_map, mechanism, mechanism_state, purview, EMD_metric, unconstrained_distributions)
        irreducibility_list[i] = irreducibility
    end # for purview

    max_irr = max(irreducibility_list...)

    if max_irr > 0
        cause_index = findall(isequal(max_irr), irreducibility_list)
        causes = purview_list[cause_index]
        cause = causes[end]
    else
        cause = [0]
    end

    return cause, max_irr

end # function get_cause

function get_system_causes(TPM, state_map, system_state, EMD_metric, unconstrained_distributions, power_set, exclude_mechanism = [], exclude_purview = [])

    n = state_map.n

    if !(isempty(exclude_mechanism))
        mechanism_list = unique([setdiff(i, exclude_mechanism) for i in 1:n])
        mechanism_list = [i for i in mechanism_list if !(isempty(i))]
    else
        mechanism_list = power_set
    end

    if !(isempty(exclude_purview))
        purview_list = unique([setdiff(i, exclude_purview) for i in power_set])
        purview_list = [i for i in purview_list if !(isempty(i))]
    else
        purview_list = power_set
    end

    n_mechanisms = length(mechanism_list)
    system_causes = [Vector() for i in 1:n_mechanisms]
    causal_strengths = zeros(n_mechanisms)

    for i in 1:n_mechanisms
        mechanism = mechanism_list[i]
        mechanism_state = system_state[mechanism]
        cause, strength = get_cause(TPM, state_map, mechanism, mechanism_state, purview_list, EMD_metric, unconstrained_distributions)
        system_causes[i] = cause
        causal_strengths[i] = strength
    end


    return system_causes, causal_strengths
end

function get_trial_streams(n_nodes, trials, animat_data, animat_TPM)

    inactive_elements = []
    for i in 3:(n_nodes+2)
        sum(animat_data[:,i]) == 0 ? push!(inactive_elements, i-2) : nothing
    end
    
    state_map = get_state_map(n_nodes)

    data = animat_data[animat_data[:,"trial"] .∈ (trials,), 1:10]

    TPM = convert_TPM(state_map, animat_TPM)

    EMD_metric = get_EMD_metric(state_map)

    power_set = collect(combinations(1:n_nodes))

    unconstrained_distributions = [[1/(2^i) for j in 1:(2^i)] for i in 1:n_nodes]

    n_trials = length(trials)
    n_timesteps = max(data[:,"timestep"]...)

    exclude_mechanism = [1,2, inactive_elements...]
    exclude_purview = [3,4, inactive_elements...]

    mechanisms = [i for i in 1:n_nodes if i ∉ exclude_mechanism]

    trial_streams = zeros(n_trials, n_timesteps, n_nodes, n_nodes)
    trial_streams_strength = zeros(n_trials, n_timesteps, n_nodes, n_nodes)

    system_state_list = []
    trial_timestep_list = []

    for row in eachrow(Matrix(data))
        
        trial = row[1]
        timestep = row[2]
        system_state = row[3:(n_nodes+2)]

        trial_index = findfirst(isequal(trial), trials)
        
        system_state_check = system_state[Not(exclude_mechanism)]

        if system_state_check in system_state_list

            copy_index = findfirst(isequal(system_state_check), system_state_list)

            copy_trial = trial_timestep_list[copy_index][1]
            copy_timestep = trial_timestep_list[copy_index][2]

            trial_streams[trial_index,timestep,:,:] = trial_streams[copy_trial,copy_timestep,:,:]
            trial_streams_strength[trial_index,timestep,:,:] = trial_streams[copy_trial,copy_timestep,:,:]
            
        else
            causes, strengths = get_system_causes(TPM, state_map, system_state, EMD_metric, unconstrained_distributions, power_set, exclude_mechanism, exclude_purview)

            for mech_index in 1:length(mechanisms)

                mech = mechanisms[mech_index]
                for cause in causes[mech_index]
                    trial_streams[trial_index,timestep,cause,mech] = 1
                    trial_streams_strength[trial_index,timestep,cause,mech] = strengths[mech_index]
                end
            end

        end
        push!(system_state_list, system_state_check)
        push!(trial_timestep_list, (trial_index,timestep))

        println((trial_index,timestep))

    end

    return trial_streams, trial_streams_strength
end 

function get_stream_perspectives_IO(trial_streams, trial_streams_strength)
    
    n_trials, n_timesteps, n_elements = size(trial_streams)

    stream_perspectives = DataFrame(trial = Int64[], persp_element = Int64[], persp_timestep = Int64[], cause_timestep = Int64[], reader_timestep = Int64[], cause = Int64[], reader = Int64[], strenght = Float64[])

    for persp_element in 1:2, trial in 1:n_trials, persp_timestep in 1:n_timesteps

        stream_check = "GO!"

        cause_list = [persp_element]
        timestep = persp_timestep

        while stream_check == "GO!"
            new_causes = []
            for cause in cause_list
                reader_list = findall(isequal(1), trial_streams[trial,timestep,cause,:])

                if reader_list != []
                    for reader in reader_list
                        strength = trial_streams_strength[trial,timestep,cause,reader]
                        row = [trial, persp_element, persp_timestep, timestep, timestep+1, cause, reader, strength]
                        push!(stream_perspectives, row)
                    end

                    push!(new_causes,setdiff(reader_list, new_causes)...)
                end
            end

            cause_list = new_causes

            new_causes == [] ? stream_check = "STOP!" : nothing
            timestep < n_timesteps ? timestep += 1 : stream_check = "STOP!"

        end

    end


    for persp_element in 3:4, trial in 1:n_trials, persp_timestep in n_timesteps:-1:1
        stream_check = "GO!"

        reader_list = [persp_element]
        timestep = persp_timestep

        while stream_check == "GO!"
            new_readers = []
            for reader in reader_list
                cause_list = findall(isequal(1), trial_streams[trial,timestep,:,reader])

                if cause_list != []
                    for cause in cause_list
                        strength = trial_streams_strength[trial,timestep,cause,reader]
                        row = [trial, persp_element, persp_timestep, timestep-1, timestep, cause, reader,strength]
                        push!(stream_perspectives, row)
                    end

                    push!(new_readers,setdiff(cause_list, new_readers)...)
                end
            end

            reader_list = new_readers
            new_readers == [] ? stream_check = "STOP!" : nothing
            timestep > 1 ? timestep -= 1 : stream_check = "STOP!"

        end

    end

    return stream_perspectives
end

