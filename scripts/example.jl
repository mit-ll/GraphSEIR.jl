# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

using Pkg
Pkg.activate("../")
using GraphSEIR
using Graphs
using MetaGraphs
using Random
using ProgressMeter
using Plots
using D3Trees
using DataStructures
pyplot()

inf_trees = SEIRInfectionTree[]

vertices_filename = "./node_attributes.json"
edges_filename = "./edge_attributes.json"
num_runs = 4
pv = 0.005
time_max = 100
scale = 10

gc = load_graph(vertices_filename, edges_filename)
num_comms = get_prop(gc, :num_communities)

#Initialize all communities with the same set of infection parameters
params = SEIRParameters([0.05, 0.0025, 0.048, 0.064]..., num_comms)

weight_params = get_prop(gc, :weight_params) * scale
set_prop!(gc, :weight_params, weight_params)

#setup compartmental analysis trackers for every run.
communities = [0] #0 is global track, add to list for other communities
trackers = [CompartmentalTracker(communities) for i in 1:num_runs]

final_nodes = [[]]
final_text = ["All Runs"]
for n in 1:num_runs
    rng = MersenneTwister(n)
    g = deepcopy(gc)
    t = expose!(g, pv, rng) #Initialize infection using built-in expose! function

    SEIR = counts(g)
    println("Run $n, initial SEIR Vals:")
    for val in SEIR
        println("$val")
    end

    e_ind = get_E_indices(g)

    for i in 1:time_max
        transition!(g=g, params=params, time_index=i, t=t, tracker=trackers[n], rng=rng)
    end
    push!(inf_trees, t)

    println("Finished run $n, Final SEIR Vals:")
    SEIR = counts(g)
    for val in SEIR
        println("$val")
    end

    function flatten!(nodes::Any, text::Any, i::Int64)
        tree = bfs_tree(inf_trees[1], i)
        q = Queue{Tuple{Int,Int,Int}}()
        # tuples represent (parent index in nodes, index in tree)
        enqueue!(q, (1, -1, i))

        while !isempty(q)
            (p, p_g, n) = dequeue!(q)
            neighb = neighbors(tree, n)
            label = string(get_prop(inf_trees[1], n, :id), "
C:", get_prop(inf_trees[1], n, :community))
            if p_g != -1
                label = string(label, "
T: ", get_prop(inf_trees[1], Edge(p_g, n), :time))
            end
            if size(neighb)[1] > 0
                label = string(label, "
(", size(neighb)[1], ")")
            end
            push!(text, label)
            push!(nodes, [])
            push!(nodes[p], size(nodes)[1])
            for v in neighb
                enqueue!(q, (size(nodes)[1], n, v))
            end
        end
    end

    nodes = [[]]
    text = ["Run $n"]
    for i in vertices(inf_trees[1])
        if parse(Int, get_prop(inf_trees[1], i, :id)) in e_ind
            flatten!(nodes, text, i)
        end
    end

    nodes = map((l) -> map((n) -> n + size(final_nodes)[1], l), nodes)

    push!(final_nodes[1], size(final_nodes)[1] + 1)
    append!(final_nodes, nodes)
    append!(final_text, text)

end

t = D3Tree(final_nodes, text=final_text)
inchrome(t)
