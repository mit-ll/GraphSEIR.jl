var documenterSearchIndex = {"docs":
[{"location":"#GraphSEIR.jl-Documentation","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"CurrentModule = GraphSEIR","category":"page"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"Disease simulation package for network diffusion SEIR model of disease transmission. In this package we use terminology networks and graphs, nodes and vertices interchangeably.","category":"page"},{"location":"#Graph-Functions","page":"GraphSEIR.jl Documentation","title":"Graph Functions","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"load_graph(vertices_filename::String, edges_filename::String)","category":"page"},{"location":"#GraphSEIR.load_graph-Tuple{String,String}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.load_graph","text":"load_graph(vertices_filename::String, edges_filename::String)\n\nReturn a MetaGraph corresponding to attributes in the corresponding vertices_filename and edges_filename.\n\n\n\n\n\n","category":"method"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"load_edgefile","category":"page"},{"location":"#GraphSEIR.load_edgefile","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.load_edgefile","text":"load_edgefile(vlist::Array{Int}, edges_filename::String; scale::Float64=1.0)\n\nReturns an edge list and weight matrix for the attributes in the corresponding edges_filename with weights scaled by option scale.\n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"swap_edges!","category":"page"},{"location":"#GraphSEIR.swap_edges!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.swap_edges!","text":"swap_edges!(g::SEIRGraph, edge_list::Dict, weight_matrix::SparseMatrixCSC)\n\nReplaces edges of graph g inplace with edges from edge_list and weights from weight_matrix.\n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"vertex_weight","category":"page"},{"location":"#GraphSEIR.vertex_weight","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.vertex_weight","text":"vertex_weight(g, v, filter_household=false)\n\nComputes weighted in-degree of a vertex v in graph g. Set filter_household flag to true to filter household edges from computation.\n\n\n\n\n\n","category":"function"},{"location":"#Disease-Functions","page":"GraphSEIR.jl Documentation","title":"Disease Functions","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"count(  g::SEIRGraph, s::SEIRStatus)\ncount(  g::SEIRGraph, c::Int, s::SEIRStatus)\ncount(  g::SEIRGraph)\ncounts(  g::SEIRGraph)","category":"page"},{"location":"#GraphSEIR.count-Tuple{MetaGraphs.MetaGraph{Int64,Float64},SEIRStatus}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.count","text":"count(  g::SEIRGraph, s::SEIRStatus)\n\nReturns the number of vertices with state s.\n\n\n\n\n\n","category":"method"},{"location":"#GraphSEIR.count-Tuple{MetaGraphs.MetaGraph{Int64,Float64},Int64,SEIRStatus}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.count","text":"count(  g::SEIRGraph, c::Int, s::SEIRStatus)\n\nReturns the number of vertices in community c that have state s.\n\n\n\n\n\n","category":"method"},{"location":"#GraphSEIR.count-Tuple{MetaGraphs.MetaGraph{Int64,Float64}}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.count","text":"count(  g::SEIRGraph)\n\nReturns the total number of vertices in the graph.\n\n\n\n\n\n","category":"method"},{"location":"#GraphSEIR.counts-Tuple{MetaGraphs.MetaGraph{Int64,Float64}}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.counts","text":"counts(  g::SEIRGraph)\n\nReturns the total number of vertices in each state S,E,I,R.\n\n\n\n\n\n","category":"method"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"counts_by_community(g::SEIRGraph, s::SEIRStatus)","category":"page"},{"location":"#GraphSEIR.counts_by_community-Tuple{MetaGraphs.MetaGraph{Int64,Float64},SEIRStatus}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.counts_by_community","text":"counts_by_community(g::SEIRGraph, s::SEIRStatus)\n\nReturns a list of the number of vertices in each community that have state s.\n\n\n\n\n\n","category":"method"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"get_indices(  g::SEIRGraph, s::SEIRStatus)\nget_indices(  g::SEIRGraph, c::Int, s::SEIRStatus)\nget_indices(  g::SEIRGraph)","category":"page"},{"location":"#GraphSEIR.get_indices-Tuple{MetaGraphs.MetaGraph{Int64,Float64},SEIRStatus}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.get_indices","text":"get_indices(  g::SEIRGraph, s::SEIRStatus)\n\nReturns an iterator of the vertices in state s.\n\n\n\n\n\n","category":"method"},{"location":"#GraphSEIR.get_indices-Tuple{MetaGraphs.MetaGraph{Int64,Float64},Int64,SEIRStatus}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.get_indices","text":"get_indices(  g::SEIRGraph, c::Int, s::SEIRStatus)\n\nReturns an iterator of all vertices in community c that have state s.\n\n\n\n\n\n","category":"method"},{"location":"#GraphSEIR.get_indices-Tuple{MetaGraphs.MetaGraph{Int64,Float64}}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.get_indices","text":"get_indices(  g::SEIRGraph)\n\nReturns iterators of the vertices in each state S,E,I,R.\n\n\n\n\n\n","category":"method"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"SEIRTransitionHistory","category":"page"},{"location":"#GraphSEIR.SEIRTransitionHistory","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.SEIRTransitionHistory","text":"SEIRTransitionHistory(SEE::T, SEI::T, EI::T, IR::T, S::T, E::T, I::T, R::T) where T <: Int\n\nStructure for tracking SEIR transition history.\n\n\n\n\n\n","category":"type"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"CompartmentalTracker(communities::Array{Int})","category":"page"},{"location":"#GraphSEIR.CompartmentalTracker-Tuple{Array{Int64,N} where N}","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.CompartmentalTracker","text":"CompartmentalTracker(communities::Array{Int})\n\nCreate dictionary for tracking compartmental and SEIR statistics for input communities.\n\n\n\n\n\n","category":"method"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"SEIRInfectionTree","category":"page"},{"location":"#GraphSEIR.SEIRInfectionTree","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.SEIRInfectionTree","text":"SEIRInfectionTree\n\nTree detailing the SEIR infection history.\n\n\n\n\n\n","category":"type"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"initialize_with!\nexpose!\nreset!","category":"page"},{"location":"#GraphSEIR.initialize_with!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.initialize_with!","text":"initialize_with!(g::SEIRGraph, p::Function, rng::AbstractRNG)\n\nInitialize statuses of vertices in graph g using density function p(g,v) that defines the probability of a vertex v being in each SEIR state. Updates are stored in infection tree t.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.expose!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.expose!","text":"expose!(g::SEIRGraph, p::AbstractFloat, rng::AbstractRNG)\n\nRandomly initialize vertices to statuse Exposed(E) with probability p.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.reset!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.reset!","text":"reset!(g::SEIRGraph, rng::AbstractRNG)\n\nSet the status of all vertices to Susceptible(S). \n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"expose_communities!\nexpose_communities","category":"page"},{"location":"#GraphSEIR.expose_communities!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.expose_communities!","text":"expose_communities!(g::SEIRGraph, p::AbstractFloat, nc::Int, rng::AbstractRNG)\n\nExpose p fraction of individuals in the entire graph, localized to nc communities.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.expose_communities","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.expose_communities","text":"expose_communities(g::SEIRGraph, v::Int, p::AbstractFloat, communities::Vector{Int})\n\nGraph-vertex function for initializing graph with exposed communities.\n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"add_new_infection!","category":"page"},{"location":"#GraphSEIR.add_new_infection!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.add_new_infection!","text":"add_new_infection!(t::SEIRInfectionTree, time_index::Int64, src::Dict, dst::Dict)\n\nAdd a new infection to the infection tree t given the src and dst vertices at time_index.\n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"transition!\nfull_transition!\nshort_circuit_transition!\nvectorized_transition!\nmdp_transition!","category":"page"},{"location":"#GraphSEIR.transition!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.transition!","text":"transition!(;g::SEIRGraph, params::Union{Nothing, SEIRParameters} = nothing, rng::AbstractRNG, tracker::Union{Nothing, Dict} = nothing, t::Union{Nothing, SEIRInfectionTree} = nothing, time_index::Int = 0)\n\nTransition graph g one time step based on the input SEIR parameters and tracking options. GlobalParameters and CommunityParameters types use the fast vectorized transition.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.full_transition!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.full_transition!","text":"full_transition!(g::SEIRGraph, params::SEIRParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict}, t::SEIRInfectionTree, time_index::Int)\n\nTransition the SEIRGraph g one time step and fully track infection tree.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.short_circuit_transition!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.short_circuit_transition!","text":"short_circuit_transition!(g::SEIRGraph, params::CommunityPairParameters, rng::AbstractRNG, tracker::Union{Nothing,Dict} = nothing)\n\nTransition the SEIRGraph g one time step using CommunityPairParameters with short circuit evaluation. Cannot track an infection tree.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.vectorized_transition!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.vectorized_transition!","text":"vectorized_transition!(g::SEIRGraph, rng::AbstractRNG)\n\nTransition the SEIRGraph g using GlobalParameters or CommunityParameters. Cannot track infection type or infection tree.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.mdp_transition!","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.mdp_transition!","text":"mdp_transition!(;g::SEIRGraph, rng::AbstractRNG)\n\nSpecialized transition function for use in EpidemicGraphMDPs planning.\n\n\n\n\n\n","category":"function"},{"location":"#Community-Ordering-Heuristics","page":"GraphSEIR.jl Documentation","title":"Community Ordering Heuristics","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"transition_likelihood\ncommunity_expected_exposure\ncommunity_expected_infection\ncommunity_expected_recovery\ncommunity_expected_exposures","category":"page"},{"location":"#GraphSEIR.transition_likelihood","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.transition_likelihood","text":"transition_likelihood(g::SEIRGraph)\n\nReturn the likelihood of all nodes transitioning to the following SEIR status given the node's current SEIR status, neighbors, and SEIR parameters.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.community_expected_exposure","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.community_expected_exposure","text":"community_expected_exposure(g::GraphSEIR, c::Int)\n\nReturn the expected number of individuals in community c that will become exposed at the next time step. (complex binomial - MLE)\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.community_expected_infection","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.community_expected_infection","text":"community_expected_infection(g::GraphSEIR, c::Int)\n\nReturn the expected number of individuals in community c that will become infected at the next time step. (binomial - MLE)\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.community_expected_recovery","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.community_expected_recovery","text":"community_expected_recovery(g::GraphSEIR, c::Int)\n\nReturn the expectation of the number of individuals in community c that will become recovered at the next time step. (binomial - MLE)\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.community_expected_exposures","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.community_expected_exposures","text":"community_expected_exposures(g::GraphSEIR)\n\nReturn the expected number of individuals in each community that will become exposed at the next time step. (complex binomial - MLE) Same function naming convention exists for infections and recoveries.\n\n\n\n\n\n","category":"function"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"exposure_subgraph\nepidemic_mobility_weight_matrix\neigenvectors\neigenscore\ncommunity_eigenscores","category":"page"},{"location":"#GraphSEIR.exposure_subgraph","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.exposure_subgraph","text":"exposure_subgraph(g::SEIRGraph)\n\nReturn the subgraph corresponding to S-E & S-I edges.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.epidemic_mobility_weight_matrix","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.epidemic_mobility_weight_matrix","text":"epidemic_mobility_weight_matrix(g::SEIRGraph)\n\nGet epidemic mobility weight matrix for corresponding SEIRGraph g.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.eigenvectors","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.eigenvectors","text":"eigenvectors(g::SEIRGraph, n::Int)\n\nReturn the n leading left (u) and right (v) eigenvectors of the epidemic-mobility weight matrix of the SEIRGraph g.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.eigenscore","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.eigenscore","text":"eigenscore(e::AbstractEdge, u::Vector, v::Vector)\n\nReturn the eigenscore u(i)*v(j) for edge e: i -> j, with left u and right v eigenvectors of the epidemic-mobility weight matrix.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.community_eigenscores","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.community_eigenscores","text":"community_eigenscores(g::SEIRGraph)\n\nReturn the sum of eigenscores for all non-household edges in each community.\n\n\n\n\n\n","category":"function"},{"location":"#SEIR-Parameter-Functions","page":"GraphSEIR.jl Documentation","title":"SEIR Parameter Functions","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"SEIRParameters\nGlobalParameters\nCommunityParameters\nCommunityPairParameters","category":"page"},{"location":"#GraphSEIR.SEIRParameters","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.SEIRParameters","text":"SEIRParameters{T}\n\nAbstract type for SEIR parameters.\n\n\n\n\n\n","category":"type"},{"location":"#GraphSEIR.GlobalParameters","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.GlobalParameters","text":"GlobalParameters(βE::T, βI::T, γ::T, λ::T)\n\nStructure for global scalar SEIR parameters.\n\n\n\n\n\n","category":"type"},{"location":"#GraphSEIR.CommunityParameters","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.CommunityParameters","text":"CommunityParams(βE::Vector{T} βI::Vector{T}, γ::T, λ::T, c::Int) where T <: Real\n\nStructure for individual community-based SEIRParameters constructor given vector β values and c communities. \n\nCan be constructed with scalar β values to set all communities to the same β values.\n\n\n\n\n\n","category":"type"},{"location":"#GraphSEIR.CommunityPairParameters","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.CommunityPairParameters","text":"CommunityPairParameters(βE::Matrix{T}, βI::Matrix{T}, γ::T, λ::T)\n\nStructure for community-pair based SEIR parameters.\n\nCan be constructed with scalar β values to set all community-pair β values.\n\n\n\n\n\n","category":"type"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"R_i\ncompute_R0\nget_betas","category":"page"},{"location":"#GraphSEIR.R_i","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.R_i","text":"R_i(g, β_e, vs = []; γ = 0.2)\n\nCompute expected number of infections generated by each vertex for specified β and γ values.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.compute_R0","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.compute_R0","text":"compute_R0(g, β_e; α=1.0, γ=0.2, filter_household = false, vweights = nothing)\n\nApproximate the reproduction rate of a disease for a given graph g, infection parameters β and γ, and top proportion of nodes by degree α. vweights can be supplied if pre-computed.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.get_betas","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.get_betas","text":"get_betas(g, R0; α=1.0, γ=0.2, λ=0.2, p = 1/5, filter_household = false, range = [-5,0])\n\nFit β_e and β_i parameters for a graph g with disease reproduction rate R0. range is the upper and lower bound as an order of magnitude for fitting. \n\n\n\n\n\n","category":"function"},{"location":"#Analysis-Helper-Functions","page":"GraphSEIR.jl Documentation","title":"Analysis Helper Functions","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"format_data","category":"page"},{"location":"#GraphSEIR.format_data","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.format_data","text":"format_data(tracker, g, pv, params, save_dir::Union{Nothing, String} = nothing)\n\nFormat simulation output from tracker into global SEIR statistics and community SEIR statistics. Saves as JSON files to save_dir if provided.\n\n\n\n\n\n","category":"function"},{"location":"#Edge-Filter-Helper-Functions","page":"GraphSEIR.jl Documentation","title":"Edge Filter Helper Functions","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"is_any_community_edge\nis_within_community_edge\nis_between_community_edge\nis_household_edge\nis_nonhousehold_edge\nis_any_community_nonhousehold_edge","category":"page"},{"location":"#GraphSEIR.is_any_community_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_any_community_edge","text":"is_any_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)\n\nReturn true if at least one vertex of edge e is in community c (e is an edge related to community c).\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.is_within_community_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_within_community_edge","text":"is_within_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)\n\nReturn true if both vertices of edge e are in community c (e is an edge within community c.)\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.is_between_community_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_between_community_edge","text":"is_any_community_edge(g::SEIRGraph, e::AbstractEdge, c::Int)\n\nReturn true if only one vertex of edge e is in community c (e is an edge between community c  and another community).\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.is_household_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_household_edge","text":"is_household_edge(g::SEIRGraph, e::AbstractEdge)\n\nReturn true if edge e is a household edge.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.is_nonhousehold_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_nonhousehold_edge","text":"is_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge)\n\nReturn true if edge e is not a household edge.\n\n\n\n\n\n","category":"function"},{"location":"#GraphSEIR.is_any_community_nonhousehold_edge","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.is_any_community_nonhousehold_edge","text":"is_any_community_nonhousehold_edge(g::SEIRGraph, e::AbstractEdge, c::Int)\n\nReturn true if at least one vertex of edge e is in community c (e is an edge between community c  and another community) and the edge is not a household edge.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"GraphSEIR.jl Documentation","title":"Index","text":"","category":"section"},{"location":"","page":"GraphSEIR.jl Documentation","title":"GraphSEIR.jl Documentation","text":"","category":"page"}]
}
