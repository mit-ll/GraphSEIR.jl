# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

using GraphSEIR
using MetaGraphs
using Graphs
using Random
using YAML
using JSON

function save_data(tracker, N, pv, params, filename)
    data = Dict("pop_size" => N, "init_inf_prop" => pv, "GraphSEIR_params" => params)
    comp_params = approximate_SEIR_params(tracker)[0]
    temp = tracker[0]
    data = merge(data, Dict("compart_params" => Dict("beta_E" => [c.βE for c in comp_params],
        "beta_I" => [c.βI for c in comp_params],
        "gamma" => [c.γ for c in comp_params],
        "lambda" => [c.λ for c in comp_params])))

    trajec = Dict("S" => [t.S for t in temp], "E" => [t.E for t in temp],
        "I" => [t.I for t in temp], "R" => [t.R for t in temp],
        "SE_E" => [t.SEE for t in temp], "SE_I" => [t.SEI for t in temp],
        "EI" => [t.EI for t in temp], "IR" => [t.IR for t in temp])
    data = merge(data, Dict(trajec))
    json_string = JSON.json(data)
    open(filename, "w") do f
        write(f, json_string)
    end
end

function expose_comms(g::SEIRGraph, v::Int64, p::AbstractFloat, communities)
    if get_prop(g, v, :community) in communities
        return ([1 - p, p, 0, 0])
    else
        return ([1, 0, 0, 0])
    end
end

function rand_community_exp(g::SEIRGraph, rng::AbstractRNG; num_comms=6, num_runs=5, pv=0.005, time_max=150, params=[0.05, 0.0025, 0.048, 0.064], track_comms=[0], file_prefix="random_community_data")
    C = get_prop(g, :num_communities)
    param_dict = Dict("beta_E" => params[1], "beta_I" => params[2], "gamma" => params[3], "lambda" => params[4])
    seirParams = SEIRParameters(params..., C)
    data = []
    communities = rand(rng, 1:C, num_comms)
    size_comms = sum(get_prop(g, :community_sizes)[communities])
    p = nv(g) * pv / size_comms
    for n in 1:num_runs
        tracker = CompartmentalTracker(track_comms)
        t = initialize_with!(g, (g, v) -> expose_comms(g, v, p, communities), rng)

        for i in 1:time_max
            transition!(g=g, params=seirParams, time_index=i, tracker=tracker, rng=rng)
        end
        filename = file_prefix * "$(pv)_$(num_comms)_$(n)"
        d = save_data(tracker, nv(g), pv, param_dict, filename)
        push!(data, d)
    end
    return (data)
end

function uniform_exp(g::SEIRGraph, rng::AbstractRNG; pv=0.005, num_runs=5, time_max=150, params=[0.05, 0.0025, 0.048, 0.064], track_comms=[0], file_prefix="uniform_data")
    param_dict = Dict("beta_E" => params[1], "beta_I" => params[2], "gamma" => params[3], "lambda" => params[4])
    C = get_prop(g, :num_communities)
    seirParams = SEIRParameters(params..., C)
    data = []
    for n in 1:num_runs
        tracker = CompartmentalTracker(track_comms)
        t = expose!(g, pv, rng)
        for i in 1:time_max
            transition!(g=g, params=seirParams, time_index=i, tracker=tracker, rng=rng)
        end
        filename = file_prefix * "$(pv)_$(n)"
        d = save_data(tracker, nv(g), pv, param_dict, filename)
        push!(data, d)
    end
    return (data)
end

function run(node_file, edge_file, experiments, expargs; scale=10, seed=1)
    gc = load_graph(node_file, edge_file)
    g = deepcopy(gc)
    weight_params = get_prop(g, :weight_params) * scale
    rng = MersenneTwister(seed)
    set_prop!(g, :weight_params, weight_params)
    for (experiment, args) in zip(experiments, expargs)
        data = experiment(g, rng; args...)
    end
end

#This part will need to be changed for each experiment. Maybe find an automated way to pass a config file of experiment arguments/experiments to be run in general.

config = YAML.load_file("config.yml")
for exp_config in config["experiments"]
    num_runs = config["globals"]["num_runs"]
    time_max = exp_config["time_max"]
    pvs = exp_config["pvs"]
    params = exp_config["params"]
    num_comms = exp_config["num_comms"]
    exprargs = []
    experiments = []
    for pv in pvs
        temp = Dict(:pv => pv, :time_max => time_max, :num_runs => num_runs, :params => params, :file_prefix => "trajectory_files/uniform_data")
        push!(exprargs, temp)
        push!(experiments, uniform_exp)
        for n in num_comms
            temp = Dict(:pv => pv, :time_max => time_max, :num_runs => num_runs, :params => params, :num_comms => n, :file_prefix => "trajectory_files/random_community_data")
            push!(exprargs, temp)
            push!(experiments, rand_community_exp)
        end
    end
    node_file = "node_attributes.json"
    edge_file = "edge_attributes.json"

    run(node_file, edge_file, experiments, exprargs)
end
