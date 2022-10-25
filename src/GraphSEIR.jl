# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

module GraphSEIR

using Revise
using Parameters
using Random
using Distributions

using StatsBase
using Statistics
using JSON: parsefile, json

using Graphs
using MetaGraphs

using LinearAlgebra
using StaticArrays
using SparseArrays
using Arpack
using Roots

const SEIRGraph = MetaGraph{Int64,Float64}
"""
    SEIRParameters{T}
Abstract type for SEIR parameters.
"""
abstract type SEIRParameters{T} end

# overload `Base.getindex` for `SEIRParameters` (single-index)
function Base.getindex(s::SEIRParameters, i::Int)

    1 <= i <= 4 || throw(BoundsError(s, i))

    if i == 1
        return s.βE
    elseif i == 2
        return s.βI
    elseif i == 3
        return s.γ
    elseif i == 4
        return s.λ
    end

end


# overload `Base.getindex` for `SEIRParameters` (multi-index)
Base.getindex(s::SEIRParameters, I::AbstractArray{Int}) = [s[i] for i in I]

include("graphs.jl")

export
    load_graph,
    load_edgefile,
    swap_edges!,
    vertex_weight

include("params.jl")
export
    GlobalParameters,
    CommunityParameters,
    CommunityPairParameters,
    R_i,
    compute_R0,
    get_betas

include("disease.jl")

export
    SEIRStatus,
    SEIRGraph,
    count,
    counts,
    count_S,
    count_E,
    count_I,
    count_R,
    count_N,
    counts_by_community,
    count_S_by_community,
    count_E_by_community,
    count_I_by_community,
    count_R_by_community,
    count_N_by_community,
    get_indices,
    get_S_indices,
    get_E_indices,
    get_I_indices,
    get_R_indices,
    SEIRTransitionHistory,
    CompartmentalTracker,
    SEIRInfectionTree,
    initialize_with!,
    expose!,
    infect!,
    reset!,
    expose_communities!,
    expose_communities,
    add_new_infection!,
    transition!,
    short_circuit_transition!,
    vectorized_transition!,
    full_transition!,
    transition_likelihood,
    community_expected_exposure,
    community_expected_exposures,
    community_expected_infection,
    community_expected_infections,
    community_expected_recovery,
    community_expected_recoveries,
    exposure_subgraph,
    epidemic_mobility_weight_matrix,
    eigenvectors,
    eigenscore,
    community_eigenscore,
    community_eigenscores,
    mdp_transition!


include("actions.jl")

export
    is_any_community_edge,
    is_within_community_edge,
    is_between_community_edge,
    is_household_edge,
    is_nonhousehold_edge,
    is_any_community_nonhousehold_edge


include("analysis.jl")

export
    CompartParams,
    approx_comp_params,
    approximate_SEIR_params,
    format_data,
    integrated_var,
    shannon_index_incidence,
    approx_dispersion,
    find_peaks

end # module
