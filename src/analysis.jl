# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

struct CompartParams <: FieldVector{4,Float64}
    βE
    βI
    γ
    λ
    R0
end

function approx_comp_params(h::SEIRTransitionHistory)
    N = sum(h[5:8])
    βE = (h.SEE * N) / (h.S * h.E)
    βI = (h.SEI * N) / (h.S * h.I)
    γ = (h.EI / h.E)
    λ = (h.IR / h.I)
    R0 = βE / γ + βI / λ
    return CompartParams(βE, βI, γ, λ, R0)
end


function approximate_SEIR_params(tracker::Dict; window=1)
    data = Dict()
    for c in keys(tracker)
        t = tracker[c]
        N = length(t)
        data[c] = [approx_comp_params(sum(t[i:i+window-1]) / window) for i in 1:N-window+1]
    end
    return (data)
end

"""
    format_data(tracker, g, pv, params, save_dir::Union{Nothing, String} = nothing)
Format simulation output from `tracker` into global SEIR statistics and community SEIR statistics. Saves as JSON files to `save_dir` if provided.
"""
function format_data(tracker, g, pv, params, save_dir::Union{Nothing,String}=nothing)
    N = nv(g)
    preamble = Dict("pop_size" => N, "init_inf_prop" => pv, "GraphSEIR_params" => params)
    data = deepcopy(preamble)
    for c in keys(tracker)
        if c == 0
            n = N
        else
            n = get_prop(g, :community_sizes)[c]
        end
        temp = tracker[c]
        d = Dict("community_size" => n,
            "S" => [t.S for t in temp], "E" => [t.E for t in temp],
            "I" => [t.I for t in temp], "R" => [t.R for t in temp],
            "SEE" => [t.SEE for t in temp], "SEI" => [t.SEI for t in temp], "EI" => [t.EI for t in temp], "IR" => [t.IR for t in temp])
        if !isnothing(save_dir)
            if c == 0
                save_path = save_dir * "/global_SEIR.json"
                save_data = merge(preamble, d)
            else
                save_path = save_dir * "/community_$(c)_SEIR.json"
                save_data = d
            end
            json_string = json(save_data)
            open(save_path, "w") do f
                write(f, json_string)
            end
        end
        data["$(c)"] = d
    end

    return data
end

function find_peaks(trajectories, window=0)
    max_time = size(trajectories, 1)
    temp = findmax(trajectories, dims=1)[2]
    peak_timeranges = []
    for t in temp
        peak_time = t[1]
        if peak_time <= window
            a = 1
        elseif peak_time + window - 1 > max_time
            a = max_time - window * 2 + 1
        else
            a = peak_time - window
        end
        b = a + window * 2 - 1
        push!(peak_timeranges, (peak_time, a, b))
    end

    return peak_timeranges
end

function shannon_index_incidence(trajectory)
    #From Samuel Scarpino's Nature paper: https://www.nature.com/articles/s41591-020-1104-0
    total_cases = sum(trajectory)
    norm_trajectory = trajectory ./ total_cases
    filter!(x -> x > 0, norm_trajectory)
    temp = -1 * sum([pij * log(pij) for pij in norm_trajectory])
    return (1 / temp)
end

function approx_dispersion(trajectories, start, finish=-1)
    if finish == -1
        finish = size(trajectories, 1)
    end
    avg_trajectory = mean(trajectories, dims=2)
    individual_sii = []
    for i in 1:size(trajectories, 2)
        sii = shannon_index_incidence(trajectories[start:finish, i])
        push!(individual_sii, sii)
    end
    avg_sii = shannon_index_incidence(avg_trajectory[start:finish])
    return (individual_sii, avg_sii)
end

function integrated_var(trajectories::Array, normalize=false, num=-1)
    traj_var = [var(t, dims=2) for t in trajectories]
    if normalize
        traj_mean = [mean(t, dims=2) for t in trajectories]
        traj_var = [v ./ m for (v, m) in zip(traj_var, traj_mean)]
    end
    if num == -1
        return [[sum(filter(!isnan, v[1:i])) for i in 1:length(v)] for v in traj_var]
    else
        return [[sum(filter(!isnan, v[1:i])) for i in 1:num] for v in traj_var]
    end
end

#=function approx_dispersion(trajectories, window = 20)
    peaks = find_peaks(trajectories, window)
    avg_trajectory = mean(trajectories, dims = 2)
    avg_peaks = find_peaks(avg_trajectory, window)
    individual_sii = []
    for i in 1:size(trajectories, 2)
        start = peaks[i][2]
        finish = peaks[i][3]
        sii = shannon_index_incidence(trajectories[start:finish,i])
        push!(individual_sii, sii)
    end
    start = avg_peaks[1][2]
    finish = avg_peaks[1][3]
    avg_sii = shannon_index_incidence(avg_trajectory[start:finish])
    indiv_peaks = [p[1] for p in peaks]
    avg_peak = avg_peaks[1][1]
    return(individual_sii, avg_sii, (indiv_peaks, avg_peak))
end=#
