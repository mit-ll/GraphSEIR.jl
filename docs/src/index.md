# GraphSEIR.jl Documentation
```@meta
CurrentModule = GraphSEIR
```
Disease simulation package for the network-SEIR model. In this package we use networks and graphs, nodes and vertices terminology interchangeably.

## Graph Functions
```@docs
load_graph(vertices_filename::String, edges_filename::String)
```

```@docs
load_edgefile
```

```@docs
swap_edges!
```

```@docs
vertex_weight
```

## Disease Functions

```@docs
count(  g::SEIRGraph, s::SEIRStatus)
count(  g::SEIRGraph, c::Int, s::SEIRStatus)
count(  g::SEIRGraph)
counts(  g::SEIRGraph)
```

```@docs
counts_by_community(g::SEIRGraph, s::SEIRStatus)
```

```@docs
get_indices(  g::SEIRGraph, s::SEIRStatus)
get_indices(  g::SEIRGraph, c::Int, s::SEIRStatus)
get_indices(  g::SEIRGraph)
```

```@docs
SEIRTransitionHistory
```

```@docs
CompartmentalTracker(communities::Array{Int})
```

```@docs
SEIRInfectionTree
```

```@docs
initialize_with!
expose!
reset!
```

```@docs
expose_communities!
expose_communities
```

```@docs
add_new_infection!
```

```@docs
transition!
full_transition!
short_circuit_transition!
vectorized_transition!
mdp_transition!
```

## Community Ordering Heuristics
```@docs
transition_likelihood
community_expected_exposure
community_expected_infection
community_expected_recovery
community_expected_exposures
```

```@docs
exposure_subgraph
epidemic_mobility_weight_matrix
eigenvectors
eigenscore
community_eigenscores
```

## SEIR Parameter Functions
```@docs
SEIRParameters
GlobalParameters
CommunityParameters
CommunityPairParameters
```

```@docs
R_i
compute_R0
get_betas
```

## Analysis Helper Functions
```@docs
format_data
```

## Edge Filter Helper Functions
```@docs
is_any_community_edge
is_within_community_edge
is_between_community_edge
is_household_edge
is_nonhousehold_edge
is_any_community_nonhousehold_edge
```

## Index

```@index
```
