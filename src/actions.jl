# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

## TODO move to edge_filters.jl

"""
    is_any_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)
Return `true` if at least one vertex of edge `e` is in community `c` (`e` is an
edge related to community `c`).
"""
is_any_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int) = get_prop(g, src(e), :community) == c || get_prop(g, dst(e), :community) == c


"""
    is_within_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)
Return `true` if both vertices of edge `e` are in community `c` (`e` is an edge
within community `c`.)
"""
is_within_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int) = get_prop(g, src(e), :community) == c && get_prop(g, dst(e), :community) == c


"""
    is_any_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)
Return `true` if only one vertex of edge `e` is in community `c` (`e` is an edge
between community `c`  and another community).
"""
is_between_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int) = get_prop(g, src(e), :community) == c ⊻ get_prop(g, dst(e), :community) == c


"""
    is_household_edge(g::SEIRGraph, e::AbstractEdge)
Return true if edge `e` is a household edge.
"""
is_household_edge(g::SEIRGraph, e::AbstractEdge) = get_prop(g, e, :type) == "household"


"""
    is_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge)
Return true if edge `e` is not a household edge.
"""
is_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge) = get_prop(g, e, :type) != "household"


"""
    is_any_community_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge, c::Int)
Return `true` if at least one vertex of edge `e` is in community `c` (`e` is an
edge between community `c`  and another community) and the edge is not a
household edge.
"""
is_any_community_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge, c::Int) = is_nonhousehold_edge(g, e) && is_any_community_edge(g, e, c)
