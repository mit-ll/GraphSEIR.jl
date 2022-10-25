# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

#]activate .
using GraphSEIR
#]activate
using Graphs
using MetaGraphs
using Random
using ProgressMeter

rng = MersenneTwister(1)

g = load_graph(vertices_filename="./node_attributes.json",
    edges_filename="./edge_attributes.json")
gc = deepcopy(g)
nvg = nv(g)
params = Array(transpose(repeat([0.02, 0.0018, 0.048, 0.064, 0], 1, get_prop(g, :num_communities))))

using Plots
pyplot()

function get_SEIR_arrays(SEIR_counts)

    t_hist = 0:length(SEIR_counts)-1
    S_hist = [SEIR[1] for SEIR in SEIR_counts]
    E_hist = [SEIR[2] for SEIR in SEIR_counts]
    I_hist = [SEIR[3] for SEIR in SEIR_counts]
    R_hist = [SEIR[4] for SEIR in SEIR_counts]

    return t_hist, S_hist, E_hist, I_hist, R_hist

end

cols = get_color_palette(:tab10, 0)

plot(size=(400, 400), box=:on)
xlabel!("Timesteps")
ylabel!("Individuals")

@showprogress for i in 1:10

    g = deepcopy(gc)
    t = GraphSEIR.expose(g, 0.35, rng)

    SEIR_counts = []
    push!(SEIR_counts, counts(g))

    for i in 1:100

        transition!(g, params, t, i, rng)
        push!(SEIR_counts, counts(g))

    end

    t_hist, S_hist, E_hist, I_hist, R_hist = get_SEIR_arrays(SEIR_counts)

    if i == 1
        plot!(t_hist, S_hist, label="S", c=cols[1], width=1, alpha=0.3)
        plot!(t_hist, E_hist, label="E", c=cols[2], width=1, alpha=0.3)
        plot!(t_hist, I_hist, label="I", c=cols[3], width=1, alpha=0.3)
        plot!(t_hist, R_hist, label="R", c=cols[4], width=1, alpha=0.3)
    else
        plot!(t_hist, S_hist, label=:none, c=cols[1], width=1, alpha=0.3)
        plot!(t_hist, E_hist, label=:none, c=cols[2], width=1, alpha=0.3)
        plot!(t_hist, I_hist, label=:none, c=cols[3], width=1, alpha=0.3)
        plot!(t_hist, R_hist, label=:none, c=cols[4], width=1, alpha=0.3)
    end

end

plot!()
