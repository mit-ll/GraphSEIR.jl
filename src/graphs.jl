# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

"""
    load_graph(vertices_filename::String, edges_filename::String)
Return a MetaGraph corresponding to attributes in the corresponding
`vertices_filename` and `edges_filename`.
"""
function load_graph(vertices_filename::String, edges_filename::String; scale::Float64=1.0)

    vertices = parsefile(vertices_filename)
    edges = parsefile(edges_filename)

    g = SEIRGraph()
    defaultweight!(g, 0)

    ## map down the node numbers and community numbers
    # # adjust "block" and "id" in the vertices
    # # adjust "source_community", "target_community", "source", and "target" in the edges
    vlist = sort(unique(vertex["id"] for vertex in vertices))
    clist = sort(unique(vertex["block"] for vertex in vertices))

    # 1-5 need vlist mappings for
    # 4 and 5 don't need clist mappings, but 1-3 do (Neela is fixing 1-3 so that they don't need mappings)

    # TODO we can also read in the the mapping json file

    for vertex in vertices

        vertex["id"] = findfirst(isequal(vertex["id"]), vlist) - 1
        vertex["block"] = findfirst(isequal(vertex["block"]), clist) - 1

    end

    for edge in edges

        edge["source"] = findfirst(isequal(edge["source"]), vlist) - 1
        edge["target"] = findfirst(isequal(edge["target"]), vlist) - 1
        # edge["source_community"] = findfirst(isequal(edge["source_community"]), clist) - 1
        # edge["target_community"] = findfirst(isequal(edge["target_community"]), clist) - 1

    end

    # number of communities
    nc = maximum((vertex["block"] for vertex in vertices)) + 1
    𝒞 = 1:nc
    𝒞map = [vertex["block"] + 1 for vertex in vertices]
    𝒞sizes = [Base.count(vertex -> vertex["block"] + 1 == c, vertices) for c in 𝒞] # TODO fix Base.count overloading

    # intialize community weight matrix and household weight
    Wc = zeros(nc, nc)
    h_idx = findfirst(e -> e["edge_type"] == "household", edges)
    wh = isnothing(h_idx) ? 1.0 : edges[h_idx]["weight"] * scale


    # construct vertices
    for vertex in vertices

        add_vertex!(g, Dict(:id => "$(vertex["id"] + 1)",
            :community => vertex["block"] + 1,
            :household => haskey(vertex, "household_id") ? vertex["household_id"] : -1,
            :status => S))

    end

    I = [edge["source"] + 1 for edge in edges]
    J = [edge["target"] + 1 for edge in edges]
    V = [Float64(edge["weight"]) for edge in edges] .* scale
    weight_mat = sparse(vcat(I, J), vcat(J, I), vcat(V, V), length(vertices), length(vertices), +)

    # construct edges
    for edge in edges

        i, j = edge["source"] + 1, edge["target"] + 1

        add_edge!(g, i, j, Dict(:type => edge["edge_type"],
            :weight => edge["weight"] * scale))

        if edge["edge_type"] != "household"

            ci, cj = get_prop(g, i, :community), get_prop(g, j, :community)

            if Wc[ci, cj] == 0.0

                Wc[ci, cj], Wc[cj, ci] = edge["weight"] * scale, edge["weight"] * scale

            end

        end

    end

    set_prop!(g, :num_communities, nc)
    set_prop!(g, :weight_params, Wc)
    set_prop!(g, :house_weight, wh)
    set_prop!(g, :vertex_community_indices, 𝒞map)
    set_prop!(g, :community_sizes, 𝒞sizes)
    set_prop!(g, :initial_adjacency, adjacency_matrix(g))
    set_prop!(g, :community_graph, Graph(Wc))
    set_prop!(g, :vlist, vlist)
    set_prop!(g, :weight_matrix, weight_mat)

    return g

end



"""
    load_edgefile(vlist::Array{Int}, edges_filename::String; scale::Float64=1.0)
Returns an edge list and weight matrix for the attributes in the corresponding `edges_filename` with weights scaled by option `scale`.
"""
function load_edgefile(vlist::Array{Int}, edges_filename::String; scale::Float64=1.0)
    # Assumes edges for a graph with N vertices. Allows for disconnected vertices.
    edges = parsefile(edges_filename)
    for edge in edges
        edge["source"] = findfirst(isequal(edge["source"]), vlist) - 1
        edge["target"] = findfirst(isequal(edge["target"]), vlist) - 1
    end

    edge_list = Dict()
    I = [edge["source"] + 1 for edge in edges]
    J = [edge["target"] + 1 for edge in edges]
    V = [Float64(edge["weight"]) * scale for edge in edges]
    N = length(vlist)

    for edge in edges
        i, j = edge["source"] + 1, edge["target"] + 1
        edge_list[(i, j)] = Dict(:type => edge["edge_type"], :weight => edge["weight"] * scale)
    end
    weight_mat = sparse(vcat(I, J), vcat(J, I), vcat(V, V), N, N, +)
    return (edge_list, weight_mat)
end


"""
    swap_edges!(g::SEIRGraph, edge_list::Dict, weight_matrix::SparseMatrixCSC)
Replaces edges of graph `g` inplace with edges from `edge_list` and weights from `weight_matrix`.
"""
function swap_edges!(g::SEIRGraph, edge_list::Dict, weight_matrix::SparseMatrixCSC)
    # Accepts output from load_edgefile and replaces edges/weight matrix for a SEIRGraph.
    edge_iter = map(x -> Graphs.SimpleEdge(x), collect(keys(edge_list)))
    g_new = SimpleGraphFromIterator(edge_iter)
    g.graph = g_new
    new_edge_list = Dict(Graphs.SimpleEdge(k) => v for (k, v) in edge_list)
    g.eprops = new_edge_list
    set_prop!(g, :weight_matrix, weight_matrix)
end


"""
    vertex_weight(g, v, filter_household=false)
Computes weighted in-degree of a vertex `v` in graph `g`. Set `filter_household` flag to true to filter household edges from computation.
"""
function vertex_weight(g, v, filter_household=false)
    N = neighbors(g, v)
    weight = 0
    for n in N
        etype = get_prop(g, v, n, :type)
        if filter_household
            if etype != "household"
                weight += get_prop(g, v, n, :weight)
            end
        else
            weight += get_prop(g, v, n, :weight)
        end
    end
    return weight
end
