# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

using Graphs
using Distributions: Poisson, Cauchy, Pareto, Levy, Multinomial, Beta
using LinearAlgebra
using StatsBase
using Random
using NPZ
using JSON
using SparseArrays
using Plots
using GraphSEIR


function get_graph_from_nodes_edgemat(num_nodes::Int64, edge_dist_mat::SparseMatrixCSC)

    g = SimpleGraph(num_nodes)

    nzentries = findall(!iszero, edge_dist_mat)

    for nze in nzentries
        edge_idxs = Tuple(nze)
        add_edge!(g, edge_idxs[1], edge_idxs[2])
    end

    return g
end


"""
    SEIRGraph(C::Vector{Int}, W::AbstractArray)
Construct a `MetaGraph` from a community membership vector `C` and edge weight
matrix `W`. The graph is constructed sequentially by adding vertices with
properties `:community => c` and `:status => S` (`SEIRStatus S`) and edges
with properties `:weight => w` and `:type => t` (t ∈ {:within, :between}).
"""
function SEIRGraph(C::Vector{Int}, W::AbstractArray)

    g = MetaGraph()

    for c in C
        add_vertex!(g, Dict(:community => c, :status => S))
    end

    for (u, v, w) in zip(findnz(W)...)
        add_edge!(g, u, v, Dict(:weight => w, :type => get_edge_type(g, u, v)))
    end

    defaultweight!(g, 0.0)

    return g

end


"""
    get_edge_type(g::MetaGraph, u::Int, v::Int)
Return the edge type based on whether the vertices in the edge are `:within` the
same `:community` or `:between` two different `:community`s.
"""
function get_edge_type(g::MetaGraph, u::Int, v::Int)
    return get_prop(g, u, :community) == get_prop(g, v, :community) ? :within : :between
end


function get_unweighted_masks(graph::AbstractGraph, max_distance::Real)
    p = floyd_warshall_shortest_paths(graph)
    dists = p.dists
    mask = map(d -> (d > max_distance), dists)
    dists[mask] .= 0
    return dists
end

function main(data_dir)

    community_sizes = npzread(data_dir * "/population_counts.npy")
    community_edge_probs = npzread(data_dir * "/edge_probabilities.npy")
    edge_weights = npzread(data_dir * "/weighted_top_edgemat.npy")
    edge_weights′ = sparse(UpperTriangular(D))

    sbm = StochasticBlockModel(community_sizes, edge_weights′)

    sbm_graph = get_graph_from_nodes_edgemat(sum(community_sizes), edge_weights′)
    #sbm_graph = SEIRGraph(sbm.nodemap, D′)

    mask = map(d -> (d != 0), D)
    D[mask] .= 0.5
    graph_attrib = GraphAttrib(sbm.nodemap, D)
    num_runs = 4
    max_distance = 1
    pv = 0.5
    time_max = 100

    inf_tracks = infTracking[]

    true_params = SEIRParameters(βE=1.997e-1, βI=1.615e-2, γ=4.763e-2, λ=5.904e-2 + 5.243e-3)
    true_Rt = true_params.βE / true_params.γ + true_params.βI / true_params.λ

    base_params = SEIRParameters(βE=0.02, βI=0.0018, γ=0.048, λ=0.064) # additional 0 for λD
    params = fill(base_params, length(community_sizes))

    for n in 1:num_runs

        println("Run number $n")
        inf_track = infTracking([], Tuple{Int64,Int64}[], [], [])
        rng = MersenneTwister(n)

        total_population = sum(community_sizes)

        population = SEIRPopulation(total_population)

        init_inf = rand(rng, Float64, length(population))

        population[init_inf.<pv/100.0] .= I

        show(population)

        for day in 1:time_max

            population = transition!(graph_attrib, population, params, inf_track, day, rng)

        end

        println("Final population stats")

        show(population)

        push!(inf_tracks, inf_track)

    end

    non_avg, avg = analyze_changes(inf_tracks)
    title = "beta_E(t) $(params[1].βE), True: $(round(true_params.βE, digits = 4))"
    plot_data((non_avg[1], avg[1]), title)
    labels = ["S" "E" "I" "R"]
    SEIR_stats = cat([hcat(inf_track.pop_changes...)[5:8, :] for inf_track in inf_tracks]..., dims=3)

    plot(SEIR_stats[:, :, 1]', label=labels, title="SEIR Run 1")
    window = 3
    non_avg, avg = analyze_changes(inf_tracks, window)
    title = "Rolling beta_E(t) $(params[1].βE), True: $(round(true_params.βE, digits = 4))"
    plot_data((non_avg[1], avg[1]), title)

    tree = gen_inf_tree(inf_tracks[1])
    R0, outdeg = approx_R0(tree)
    max_degree = maximum(outdeg)
    histogram(outdeg, bins=0:max_degree,
        title="Hist. of secondary infections, max: $(max_degree)",
        label="R0 = $(round(R0,digits = 3))")
end
