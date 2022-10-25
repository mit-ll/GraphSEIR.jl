# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

"""
    SEIRStatus
Enumeration of SEIR statuses.
"""
@enum SEIRStatus S = 0 E = 1 I = 2 R = 3

""" 
    count(  g::SEIRGraph, s::SEIRStatus)
Returns the number of vertices with state `s`.
"""
count(g::SEIRGraph, s::SEIRStatus) = length(get_indices(g, s))

"""
    count(  g::SEIRGraph)
Returns the total number of vertices in the graph.
"""
count(g::SEIRGraph) = sum(counts(g))

"""
    counts(  g::SEIRGraph)
Returns the total number of vertices in each state S,E,I,R.
"""
counts(g::SEIRGraph) = count_S(g), count_E(g), count_I(g), count_R(g)
count_S(g::SEIRGraph) = count(g, S)
count_E(g::SEIRGraph) = count(g, E)
count_I(g::SEIRGraph) = count(g, I)
count_R(g::SEIRGraph) = count(g, R)
count_N(g::SEIRGraph) = count(g)


# methods for getting SEIR indices (individuals)
"""
    get_indices(  g::SEIRGraph, s::SEIRStatus)
Returns an iterator of the vertices in state `s`.
"""
get_indices(g::SEIRGraph, s::SEIRStatus) = collect(filter_vertices(g, :status, s))

"""
    get_indices(  g::SEIRGraph)
Returns iterators of the vertices in each state S,E,I,R.
"""
get_indices(g::SEIRGraph) = get_S_indices(g), get_E_indices(g), get_I_indices(g), get_R_indices(g)
get_S_indices(g::SEIRGraph) = get_indices(g, S)
get_E_indices(g::SEIRGraph) = get_indices(g, E)
get_I_indices(g::SEIRGraph) = get_indices(g, I)
get_R_indices(g::SEIRGraph) = get_indices(g, R)


# methods for getting SEIR counts within community c
""" 
    count(  g::SEIRGraph, c::Int, s::SEIRStatus)
Returns the number of vertices in community `c` that have state `s`.
"""
count(g::SEIRGraph, c::Int, s::SEIRStatus) = length(get_indices(g, c, s))
count(g::SEIRGraph, c::Int) = sum(counts(g, c))
counts(g::SEIRGraph, c::Int) = count_S(g, c), count_E(g, c), count_I(g, c), count_R(g, c)
count_S(g::SEIRGraph, c::Int) = count(g, c, S)
count_E(g::SEIRGraph, c::Int) = count(g, c, E)
count_I(g::SEIRGraph, c::Int) = count(g, c, I)
count_R(g::SEIRGraph, c::Int) = count(g, c, R)
count_N(g::SEIRGraph, c::Int) = count(g, c)


# fast count method by community
"""
    counts_by_community(g::SEIRGraph, s::SEIRStatus)
Returns a list of the number of vertices in each community that have state `s`.
"""
function counts_by_community(g::SEIRGraph, s::SEIRStatus)

    vertex_community_indices = get_prop(g, :vertex_community_indices)
    status_counts = zeros(Int, get_prop(g, :num_communities))

    for v in vertices(g)

        get_prop(g, v, :status) == s ? status_counts[vertex_community_indices[v]] += 1 : nothing

    end

    return status_counts

end


count_S_by_community(g::SEIRGraph) = counts_by_community(g, S)
count_E_by_community(g::SEIRGraph) = counts_by_community(g, E)
count_I_by_community(g::SEIRGraph) = counts_by_community(g, I)
count_R_by_community(g::SEIRGraph) = counts_by_community(g, R)
count_N_by_community(g::SEIRGraph) = get_prop(g, :community_sizes)


# methods for getting SEIR indices (individuals) within community c
"""
    get_indices(  g::SEIRGraph, c::Int, s::SEIRStatus)
Returns an iterator of all vertices in community `c` that have state `s`.
"""
get_indices(g::SEIRGraph, c::Int, s::SEIRStatus) = collect(filter_vertices(g, (g, v) -> get_prop(g, v, :community) == c && get_prop(g, v, :status) == s))
get_indices(g::SEIRGraph, c::Int) = get_S_indices(g, c), get_E_indices(g, c), get_I_indices(g, c), get_R_indices(g, c)
get_S_indices(g::SEIRGraph, c::Int) = get_indices(g, c, S)
get_E_indices(g::SEIRGraph, c::Int) = get_indices(g, c, E)
get_I_indices(g::SEIRGraph, c::Int) = get_indices(g, c, I)
get_R_indices(g::SEIRGraph, c::Int) = get_indices(g, c, R)


"""
    SEIRTransitionHistory(SEE::T, SEI::T, EI::T, IR::T, S::T, E::T, I::T, R::T) where T <: Int
Structure for tracking SEIR transition history.
"""
@with_kw mutable struct SEIRTransitionHistory{T} <: FieldVector{8,T}
    SEE::T = 0
    SEI::T = 0
    EI::T = 0
    IR::T = 0
    S::T = 0
    E::T = 0
    I::T = 0
    R::T = 0
end

"""
    CompartmentalTracker(communities::Array{Int})
Create dictionary for tracking compartmental and SEIR statistics for input communities.
"""
CompartmentalTracker(communities::Array{Int}) = Dict(c => [SEIRTransitionHistory{Float64}()] for c in communities)

"""
    SEIRInfectionTree
Tree detailing the SEIR infection history.
"""
const SEIRInfectionTree = MetaDiGraph{Int64,Float64}


"""
    initialize_with!(g::SEIRGraph, p::Function, rng::AbstractRNG)
Initialize statuses of vertices in graph `g` using density function `p(g,v)`
that defines the probability of a vertex v being in each SEIR state. Updates are stored in infection tree `t`.
"""
function initialize_with!(g::SEIRGraph, p::Function, rng::AbstractRNG)

    t = SEIRInfectionTree()

    for v in vertices(g)

        s = sample(rng, [instances(SEIRStatus)...], Weights(p(g, v)))
        set_prop!(g, v, :status, s)

        if s != S
            new_infection = copy(props(g, v))
            delete!(new_infection, :status)
            add_vertex!(t, new_infection)
        end

    end

    set_indexing_prop!(t, :id)
    return t
end

"""
    expose!(g::SEIRGraph, p::AbstractFloat, rng::AbstractRNG)
Randomly initialize vertices to statuse Exposed(E) with probability `p`.
"""
expose!(g::SEIRGraph, p::AbstractFloat, rng::AbstractRNG) = initialize_with!(g, (g, v) -> [(1 - p), p, 0, 0], rng)
infect!(g::SEIRGraph, p::AbstractFloat, rng::AbstractRNG) = initialize_with!(g, (g, v) -> [(1 - p), 0, p, 0], rng)

"""
    reset!(g::SEIRGraph, rng::AbstractRNG)
Set the status of all vertices to Susceptible(S). 
"""
reset!(g::SEIRGraph, rng::AbstractRNG) = initialize_with!(g, (g, v) -> [1, 0, 0, 0], rng)


"""
    expose_communities!(g::SEIRGraph, p::AbstractFloat, nc::Int, rng::AbstractRNG)
Expose `p` fraction of individuals in the entire graph, localized to `nc`
communities.
"""
function expose_communities!(g::SEIRGraph, p::AbstractFloat, nc::Int, rng::AbstractRNG)

    cids = sample(rng, 1:get_prop(g, :num_communities), nc, replace=false)
    nv_in_cids = sum(get_prop(g, :community_sizes)[cids])

    p_in_cids = nv(g) * p / nv_in_cids

    t = initialize_with!(g, (g, v) -> expose_communities(g, v, p_in_cids, cids), rng)

    return t

end


"""
    expose_communities(g::SEIRGraph, v::Int, p::AbstractFloat, communities::Vector{Int})
Graph-vertex function for initializing graph with exposed communities.
"""
function expose_communities(g::SEIRGraph, v::Int, p::AbstractFloat, communities::Vector{Int})

    return get_prop(g, v, :community) ∈ communities ? [1 - p, p, 0, 0] : [1, 0, 0, 0]

end


"""
    add_new_infection!(t::SEIRInfectionTree, time_index::Int64, src::Dict, dst::Dict)
Add a new infection to the infection tree `t` given the `src` and `dst` vertices
at `time_index`.
"""
function add_new_infection!(g::SEIRGraph, src_id::String, dst::Dict, time_index::Int64, t::SEIRInfectionTree)

    dst_id = pop!(dst, :id)
    if isempty(filter_vertices(t, :id, src_id))
        println("Missing infection source $src_id, did not add to tree")
    else

        if isempty(filter_vertices(t, :id, dst_id))
            add_vertex!(t, :id, dst_id)
            delete!(dst, :status)

            set_props!(t, t[dst_id, :id], dst)
        end
        add_edge!(t, t[src_id, :id], t[dst_id, :id], :time, time_index)
    end
end


"""
    transition!(;g::SEIRGraph, params::Union{Nothing, SEIRParameters} = nothing, rng::AbstractRNG, tracker::Union{Nothing, Dict} = nothing, t::Union{Nothing, SEIRInfectionTree} = nothing, time_index::Int = 0)
Transition graph `g` one time step based on the input SEIR parameters and tracking options.
`GlobalParameters` and `CommunityParameters` types use the fast vectorized transition.
"""
function transition!(; g::SEIRGraph, params::Union{Nothing,SEIRParameters}=nothing, rng::AbstractRNG, tracker::Union{Nothing,Dict}=nothing, tree::Union{Nothing,SEIRInfectionTree}=nothing, time_index::Int=0)
    if isnothing(params)
        params = get_prop(g, :SEIR_parameters)
    end
    if !isnothing(tree)
        full_transition!(g, params, rng, tracker, tree, time_index)
    elseif isa(params, GlobalParameters) || isa(params, CommunityParameters)
        vectorized_transition!(g, params, rng, tracker)
    elseif isa(params, CommunityPairParameters)
        short_circuit_transition!(g, params, rng, tracker)
    end
end


"""
    full_transition!(g::SEIRGraph, params::SEIRParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int) 
Transition the `SEIRGraph` `g` one time step and fully track infection tree.
"""
function full_transition!(g::SEIRGraph, params::SEIRParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int)
    if isa(params, GlobalParameters)
        _full_transition_global!(g, params, rng, tracker, t, time_index)
    elseif isa(params, CommunityParameters)
        _full_transition_community!(g, params, rng, tracker, t, time_index)
    elseif isa(params, CommunityPairParameters)
        _full_transition_communitypair!(g, params, rng, tracker, t, time_index)
    end
end

"""
    short_circuit_transition!(g::SEIRGraph, params::CommunityPairParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict} = nothing)
Transition the `SEIRGraph` `g` one time step using CommunityPairParameters with short circuit evaluation. Cannot track an infection tree.
"""
function short_circuit_transition!(g::SEIRGraph, params::CommunityPairParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}=nothing)

    if !isnothing(tracker)
        temp = Dict()
        global_track = false
        for c in keys(tracker)
            if c == 0
                global_track = true
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g)...])
            else
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g, c)...])
            end
        end
    end
    if isnothing(params)
        params = get_prop(g, :SEIR_parameters)
    end

    βE = params[1]
    βI = params[2]
    γ = params[3]
    λ = params[4]

    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)
    for j in E_inds
        cj = get_prop(g, j, :community)
        for i in neighbors(g, j)
            if get_prop(g, i, :status) == S
                ci = get_prop(g, i, :community)
                w = get_prop(g, i, j, :weight)
                p = 1 - exp(-βE[ci, cj] * w)
                if rand(rng) < p
                    set_prop!(g, i, :status, E)
                    if !isnothing(tracker)
                        if global_track
                            temp[0][1] += 1
                        end
                        if haskey(temp, ci)
                            temp[ci][1] += 1
                        end
                    end
                end
            end
        end
        if rand(rng) < γ
            set_prop!(g, j, :status, I)
            if !isnothing(tracker)
                if global_track
                    temp[0][3] += 1
                end
                if haskey(temp, cj)
                    temp[cj][3] += 1
                end
            end
        end
    end
    for j in I_inds
        cj = get_prop(g, j, :community)
        for i in neighbors(g, j)
            if get_prop(g, i, :status) == S
                ci = get_prop(g, i, :community)
                w = get_prop(g, i, j, :weight)
                p = 1 - exp(-βI[ci, cj] * w)
                if rand(rng) < p
                    set_prop!(g, i, :status, E)
                    if !isnothing(tracker)
                        if global_track
                            temp[0][2] += 1
                        end
                        if haskey(temp, ci)
                            temp[ci][2] += 1
                        end
                    end
                end
            end
        end
        if rand(rng) < λ
            set_prop!(g, j, :status, R)
            if !isnothing(tracker)
                if global_track
                    temp[0][4] += 1
                end
                if haskey(temp, cj)
                    temp[cj][4] += 1
                end
            end
        end
    end
    if !isnothing(tracker)
        for c in keys(tracker)
            push!(tracker[c], temp[c])
        end
    end
    return nothing
end

"""
    vectorized_transition!(g::SEIRGraph, rng::AbstractRNG)
Transition the `SEIRGraph` `g` using `GlobalParameters` or `CommunityParameters`. Cannot track infection type or infection tree.
"""
function vectorized_transition!(g::SEIRGraph, params::Union{GlobalParameters,CommunityParameters}, rng::AbstractRNG, tracker::Union{Nothing,Dict}=nothing)

    if !has_prop(g, :weight_matrix)
        weight_matrix = sparse(Graphs.weights(g))
        set_prop!(g, :weight_matrix, weight_matrix)
    end

    if !isnothing(tracker)
        temp = Dict()
        global_track = false
        for c in keys(tracker)
            if c == 0
                global_track = true
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g)...])
            else
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g, c)...])
            end
        end
    end

    if isnothing(params)
        params = get_prop(g, :SEIR_parameters)
    end
    if isa(params, CommunityParameters)
        cmap = get_prop(g, :vertex_community_indices)
        βE = [params[1][cmap[v]] for v in 1:nv(g)]
        βI = [params[2][cmap[v]] for v in 1:nv(g)]
    else
        βE = params[1]
        βI = params[2]
    end
    γ = params[3]
    λ = params[4]

    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)

    S_indicator = zeros(nv(g))
    S_indicator[S_inds] .= 1.0

    E_indicator = zeros(nv(g))
    E_indicator[E_inds] .= 1.0

    I_indicator = zeros(nv(g))
    I_indicator[I_inds] .= 1.0

    W = get_prop(g, :weight_matrix)
    pressure = W * (βE .* E_indicator + βI .* I_indicator)

    prob = (1 .- exp.(-pressure)) .* S_indicator

    for i in E_inds
        prob[i] = γ
    end
    for i in I_inds
        prob[i] = λ
    end

    for i in S_inds
        if rand(rng) < prob[i]
            set_prop!(g, i, :status, E)
            ci = get_prop(g, i, :community)
            if !isnothing(tracker)
                if global_track
                    temp[0][1] += 1
                end
                if haskey(temp, ci)
                    temp[ci][1] += 1
                end
            end
        end
    end

    for i in E_inds
        if rand(rng) < γ
            set_prop!(g, i, :status, I)
            ci = get_prop(g, i, :community)
            if !isnothing(tracker)
                if global_track
                    temp[0][3] += 1
                end
                if haskey(temp, ci)
                    temp[ci][3] += 1
                end
            end
        end
    end

    for i in I_inds
        if rand(rng) < λ
            set_prop!(g, i, :status, R)
            ci = get_prop(g, i, :community)
            if !isnothing(tracker)
                if global_track
                    temp[0][4] += 1
                end
                if haskey(temp, ci)
                    temp[ci][4] += 1
                end
            end
        end
    end
    if !isnothing(tracker)
        for c in keys(tracker)
            push!(tracker[c], temp[c])
        end
    end
end


function _full_transition_global!(g::SEIRGraph, params::GlobalParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int)

    if !isnothing(tracker)
        temp = Dict()
        global_track = false
        for c in keys(tracker)
            if c == 0
                global_track = true
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g)...])
            else
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g, c)...])
            end
        end
    end

    βE = params[1]
    βI = params[2]
    γ = params[3]
    λ = params[4]

    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)
    for i in S_inds
        first_inf = -1
        inf_type = 0
        inf_roll = 1.1
        ci = get_prop(g, i, :community)
        for j in neighbors(g, i)
            state_j = get_prop(g, j, :status)
            cj = get_prop(g, j, :community)
            w = get_prop(g, i, j, :weight)
            if state_j == E
                p = 1 - exp(-βE * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 1
                    end
                end
            elseif state_j == I
                p = 1 - exp(-βI * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 2
                    end
                end
            end
        end
        if first_inf != -1
            source = get_prop(g, first_inf, :id)
            target = copy(props(g, i))
            add_new_infection!(g, source, target, time_index, t)
        end
        if !isnothing(tracker) && first_inf != -1
            if global_track
                temp[0][inf_type] += 1
            end
            if haskey(temp, ci)
                temp[ci][inf_type] += 1
            end
        end
    end
    for i in E_inds
        if rand(rng) < γ
            set_prop!(g, i, :status, I)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][3] += 1
                end
                if haskey(temp, ci)
                    temp[ci][3] += 1
                end
            end
        end
    end
    for i in I_inds
        if rand(rng) < λ
            set_prop!(g, i, :status, R)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][4] += 1
                end
                if haskey(temp, ci)
                    temp[ci][4] += 1
                end
            end
        end
    end
    if !isnothing(tracker)
        for c in keys(tracker)
            push!(tracker[c], temp[c])
        end
    end
end


function _full_transition_community!(g::SEIRGraph, params::CommunityParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int)

    if !isnothing(tracker)
        temp = Dict()
        global_track = false
        for c in keys(tracker)
            if c == 0
                global_track = true
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g)...])
            else
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g, c)...])
            end
        end
    end

    βE = params[1]
    βI = params[2]
    γ = params[3]
    λ = params[4]

    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)
    for i in S_inds
        first_inf = -1
        inf_type = 0
        inf_roll = 1.1
        ci = get_prop(g, i, :community)
        for j in neighbors(g, i)
            state_j = get_prop(g, j, :status)
            cj = get_prop(g, j, :community)
            w = get_prop(g, i, j, :weight)
            if state_j == E
                p = 1 - exp(-βE[cj] * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 1
                    end
                end
            elseif state_j == I
                p = 1 - exp(-βI[cj] * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 2
                    end
                end
            end
        end
        if first_inf != -1
            source = get_prop(g, first_inf, :id)
            target = copy(props(g, i))
            add_new_infection!(g, source, target, time_index, t)
        end
        if !isnothing(tracker) && first_inf != -1
            if global_track
                temp[0][inf_type] += 1
            end
            if haskey(temp, ci)
                temp[ci][inf_type] += 1
            end
        end
    end
    for i in E_inds
        if rand(rng) < γ
            set_prop!(g, i, :status, I)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][3] += 1
                end
                if haskey(temp, ci)
                    temp[ci][3] += 1
                end
            end
        end
    end
    for i in I_inds
        if rand(rng) < λ
            set_prop!(g, i, :status, R)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][4] += 1
                end
                if haskey(temp, ci)
                    temp[ci][4] += 1
                end
            end
        end
    end
    if !isnothing(tracker)
        for c in keys(tracker)
            push!(tracker[c], temp[c])
        end
    end
end

function _full_transition_communitypair!(g::SEIRGraph, params::CommunityPairParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int)

    if !isnothing(tracker)
        temp = Dict()
        global_track = false
        for c in keys(tracker)
            if c == 0
                global_track = true
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g)...])
            else
                temp[c] = SEIRTransitionHistory([zeros(4)..., counts(g, c)...])
            end
        end
    end

    βE = params[1]
    βI = params[2]
    γ = params[3]
    λ = params[4]

    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)
    for i in S_inds
        first_inf = -1
        inf_type = 0
        inf_roll = 1.1
        ci = get_prop(g, i, :community)
        for j in neighbors(g, i)
            state_j = get_prop(g, j, :status)
            cj = get_prop(g, j, :community)
            w = get_prop(g, i, j, :weight)
            if state_j == E
                p = 1 - exp(-βE[ci, cj] * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 1
                    end
                end
            elseif state_j == I
                p = 1 - exp(-βI[ci, cj] * w)
                r = rand(rng)
                if r < p
                    set_prop!(g, i, :status, E)
                    if r < inf_roll
                        inf_roll = r
                        first_inf = j
                        inf_type = 2
                    end
                end
            end
        end
        if first_inf != -1
            source = get_prop(g, first_inf, :id)
            target = copy(props(g, i))
            add_new_infection!(g, source, target, time_index, t)
        end
        if !isnothing(tracker) && first_inf != -1
            if global_track
                temp[0][inf_type] += 1
            end
            if haskey(temp, ci)
                temp[ci][inf_type] += 1
            end
        end
    end
    for i in E_inds
        if rand(rng) < γ
            set_prop!(g, i, :status, I)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][3] += 1
                end
                if haskey(temp, ci)
                    temp[ci][3] += 1
                end
            end
        end
    end
    for i in I_inds
        if rand(rng) < λ
            set_prop!(g, i, :status, R)
            if !isnothing(tracker)
                ci = get_prop(g, i, :community)
                if global_track
                    temp[0][4] += 1
                end
                if haskey(temp, ci)
                    temp[ci][4] += 1
                end
            end
        end
    end
    if !isnothing(tracker)
        for c in keys(tracker)
            push!(tracker[c], temp[c])
        end
    end
end

"""
    mdp_transition!(;g::SEIRGraph, rng::AbstractRNG)
Specialized transition function for use in EpidemicGraphMDPs planning.
"""
function mdp_transition!(; g::SEIRGraph, rng::AbstractRNG)
    params = get_prop(g, :SEIR_parameters)
    cmap = get_prop(g, :vertex_community_indices)
    βE = [params[1][cmap[v]] for v in 1:nv(g)]
    βI = [params[2][cmap[v]] for v in 1:nv(g)]
    γ = params[3]
    λ = params[4]

    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)

    S_indicator = zeros(nv(g))
    S_indicator[S_inds] .= 1.0

    E_indicator = zeros(nv(g))
    E_indicator[E_inds] .= 1.0

    I_indicator = zeros(nv(g))
    I_indicator[I_inds] .= 1.0

    W = get_prop(g, :weight_matrix)
    pressure = W * (βE .* E_indicator + βI .* I_indicator)

    prob = (1 .- exp.(-pressure)) .* S_indicator

    for i in E_inds
        prob[i] = γ
    end
    for i in I_inds
        prob[i] = λ
    end

    for i in S_inds
        if rand(rng) < prob[i]
            set_prop!(g, i, :status, E)
            ci = get_prop(g, i, :community)
        end
    end

    for i in E_inds
        if rand(rng) < γ
            set_prop!(g, i, :status, I)
            ci = get_prop(g, i, :community)
        end
    end

    for i in I_inds
        if rand(rng) < λ
            set_prop!(g, i, :status, R)
            ci = get_prop(g, i, :community)
        end
    end
end

"""
    transition_likelihood(g::SEIRGraph)
Return the likelihood of all nodes transitioning to the following SEIR status
given the node's current SEIR status, neighbors, and SEIR parameters.
"""
function transition_likelihood(g::SEIRGraph, params::Union{Nothing,GlobalParameters,CommunityParameters,CommunityPairParameters}=nothing)
    if isnothing(params)
        params = get_prop(g, :SEIR_parameters)
    end
    if isa(params, GlobalParameters) || isa(params, CommunityParameters)
        return _transition_likelihood_vectorized(g, params)
    elseif isa(params, CommunityPairParameters)
        return [_transition_likelihood_pair(g, params, v) for v in vertices(g)]
    end
end

function _transition_likelihood_vectorized(g::SEIRGraph, params::Union{GlobalParameters,CommunityParameters})
    if !has_prop(g, :weight_matrix)
        weight_matrix = sparse(Graphs.weights(g))
        set_prop!(g, :weight_matrix, weight_matrix)
    end
    if isa(params, CommunityParameters)
        cmap = get_prop(g, :vertex_community_indices)
        βE = [params[1][cmap[v]] for v in 1:nv(g)]
        βI = [params[2][cmap[v]] for v in 1:nv(g)]
    else
        βE = params[1]
        βI = params[2]
    end
    γ = params[3]
    λ = params[4]

    W = get_prop(g, :weight_matrix)
    S_inds = get_S_indices(g)
    E_inds = get_E_indices(g)
    I_inds = get_I_indices(g)

    S_indicator = zeros(nv(g))
    S_indicator[S_inds] .= 1.0

    E_indicator = zeros(nv(g))
    E_indicator[E_inds] .= 1.0

    I_indicator = zeros(nv(g))
    I_indicator[I_inds] .= 1.0

    pressure = W * (βE .* E_indicator + βI .* I_indicator)

    prob = (1 .- exp.(-pressure)) .* S_indicator

    for i in E_inds
        prob[i] = γ
    end
    for i in I_inds
        prob[i] = λ
    end
    return prob
end


function _transition_likelihood_pair(g::SEIRGraph, params::CommunityPairParameters, i::Int)
    si = get_prop(g, i, :status)
    if si == S
        pressure = 0
        ci = get_prop(g, i, :community)
        for j in neighbors(g, i)
            sj = get_prop(g, j, :status)
            if j == E
                cj = get_prop(g, j, :community)
                pressure += params[1][ci, cj] * get_prop(g, i, j, :weight)
            elseif j == I
                pressure += params[2][ci, cj] * get_prop(g, i, j, :weight)
            end
        end
        return 1 - exp(-pressure)
    elseif si == E
        return params[3]
    elseif si == I
        return params[4]
    else
        return 0
    end
end


"""
    community_expected_exposure(g::GraphSEIR, c::Int)
Return the expected number of individuals in community `c` that will
become exposed at the next time step. (complex binomial - MLE)
"""
function community_expected_exposure(g::SEIRGraph, c::Int)
    if !isempty(get_S_indices(g, c))
        prob = transition_likelihood(g)
        return sum(prob[get_S_indices(g, c)])
    else
        return 0.0
    end
end


"""
    community_expected_exposures(g::GraphSEIR)
Return the expected number of individuals in each community that will
become exposed at the next time step. (complex binomial - MLE)
Same function naming convention exists for infections and recoveries.
"""
function community_expected_exposures(g::SEIRGraph)
    comm_exp_exposures = zeros(get_prop(g, :num_communities))
    cmap = get_prop(g, :vertex_community_indices)
    prob = transition_likelihood(g)
    for v in get_S_indices(g)
        cv = cmap[v]
        comm_exp_exposures[cv] += prob[v]
    end

    return comm_exp_exposures

end


"""
    community_expected_infection(g::GraphSEIR, c::Int)
Return the expected number of individuals in community `c` that will
become infected at the next time step. (binomial - MLE)
"""
function community_expected_infection(g::SEIRGraph, c::Int)
    if !isempty(get_S_indices(g, c))
        prob = transition_likelihood(g)
        return sum(prob[get_E_indices(g, c)])
    else
        return 0.0
    end
end



"""
    community_expected_infections(g::GraphSEIR)
Return the expectation of the number of individuals in each community that will
become infected at the next time step. (binomial - MLE)
"""
function community_expected_infections(g::SEIRGraph)
    comm_exp_exposures = zeros(get_prop(g, :num_communities))
    cmap = get_prop(g, :vertex_community_indices)
    prob = transition_likelihood(g)
    for v in get_E_indices(g)
        cv = cmap[v]
        comm_exp_exposures[cv] += prob[v]
    end

    return comm_exp_exposures

end


"""
    community_expected_recovery(g::GraphSEIR, c::Int)
Return the expectation of the number of individuals in community `c` that will
become recovered at the next time step. (binomial - MLE)
"""
function community_expected_recovery(g::SEIRGraph, c::Int)
    if !isempty(get_S_indices(g, c))
        prob = transition_likelihood(g)
        return sum(prob[get_I_indices(g, c)])
    else
        return 0.0
    end
end


"""
    community_expected_recoveries(g::GraphSEIR)
Return the expectation of the number of individuals in each community that will
become recovered at the next time step. (binomial - MLE)
"""
function community_expected_recoveries(g::SEIRGraph)
    comm_exp_exposures = zeros(get_prop(g, :num_communities))
    cmap = get_prop(g, :vertex_community_indices)
    prob = transition_likelihood(g)
    for v in get_I_indices(g)
        cv = cmap[v]
        comm_exp_exposures[cv] += prob[v]
    end

    return comm_exp_exposures

end


"""
    exposure_subgraph(g::SEIRGraph)
Return the subgraph corresponding to S-E & S-I edges.
"""
function exposure_subgraph(g::SEIRGraph)

    return induced_subgraph(g, filter_edges(g, (g, e) ->
        (get_prop(g, src(e), :status) == S && get_prop(g, dst(e), :status) ∈ (E, I)) ||
            (get_prop(g, dst(e), :status) == S && get_prop(g, src(e), :status) ∈ (E, I))))

end


"""
    epidemic_mobility_weight_matrix(g::SEIRGraph)
Get epidemic mobility weight matrix for corresponding `SEIRGraph` `g`.
"""
function epidemic_mobility_weight_matrix(g::SEIRGraph) # we could use only mobility weight matrix if we were doing SIR, but SEIR has two β parameters, which means S-E, S-I pathways are different contributions to epidemic

    # most of time cost comes from converting a metaweights into a matrix (sparse or dense)
    W = sparse(Graphs.weights(g)) # mobility weights

    parameters = get_prop(g, :SEIR_parameters)

    i_list, j_list, _ = findnz(W)

    # rescale mobility weight matrix by
    for (i, j) in zip(i_list, j_list)

        si, sj = get_prop(g, i, :status), get_prop(g, j, :status)

        if si == S && sj ∈ (E, I)
            ci, cj = get_prop(g, i, :community), get_prop(g, j, :community)
            W[i, j] *= parameters[Int(sj)][cj]
        elseif sj == S && si ∈ (E, I)
            ci, cj = get_prop(g, i, :community), get_prop(g, j, :community)
            W[i, j] *= parameters[Int(si)][ci]
        else
            W[i, j] *= 0
        end

    end

    # @assert issymmetric(W)
    return dropzeros(W)

end


"""
    eigenvectors(g::SEIRGraph, n::Int)
Return the `n` leading left (u) and right (v) eigenvectors of the
epidemic-mobility weight matrix of the `SEIRGraph` `g`.
"""
function eigenvectors(g::SEIRGraph, n::Int)

    W = epidemic_mobility_weight_matrix(g)

    if size(W) != (0, 0)
        svd = svds(W, nsv=n)[1] # since W is symmetric, SVD U and V give left and right eigenvectors
        return (u=svd.U, v=svd.V)
    else
        return (u=[0.0], v=[0.0])
    end

end


"""
    eigenscore(e::AbstractEdge, u::Vector, v::Vector)
Return the eigenscore u(i)*v(j) for edge e: i -> j, with left `u` and right
`v` eigenvectors of the epidemic-mobility weight matrix.
"""
eigenscore(e::AbstractEdge, u, v) = u[src(e)] * v[dst(e)]


"""
    community_eigenscore(g::SEIRGraph, c::Int, u::Vector, v::Vector)
Return the sum of eigenscores for all non-household edges in community `c`.
"""
community_eigenscore(g::SEIRGraph, c::Int, u, v) = sum(eigenscore(e, u, v) for e in filter_edges(g, (g, e) -> is_any_community_nonhousehold_edge(g, e, c)))


"""
    community_eigenscores(g::SEIRGraph)
Return the sum of eigenscores for all non-household edges in each community.
"""
function community_eigenscores(g::SEIRGraph)

    sg, vmap = exposure_subgraph(g)
    u, v = eigenvectors(sg, 1)

    all_community_eigenscores = zeros(get_prop(g, :num_communities))

    for e in filter_edges(sg, (g, e) -> is_nonhousehold_edge(g, e))

        ci = get_prop(sg, src(e), :community)
        cj = get_prop(sg, dst(e), :community)

        es = eigenscore(e, u, v)

        all_community_eigenscores[ci] += es
        all_community_eigenscores[cj] += es

    end

    return all_community_eigenscores

end
