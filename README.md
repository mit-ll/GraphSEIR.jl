# GraphSEIR.jl

GraphSEIR is a network-based disease simulation model. This model takes as input a graph or network instance that captures contact patterns in a population. The model then uses an agent-based model called SEIR (Susceptible, Exposed, Infected, Recovered) that describes how disease spreads across this graph topology to produce simulated infection statistics.

This is a Julia package and requires a local installation of [Julia](https://julialang.org). Documentation for this tool can be found [here](https://mit-ll.github.io/GraphSEIR.jl/build/).

## Global Packages
GraphSEIR makes use of "Plots" and "PyPlot" - please add those to your default Julia environment via Pkg before importing GraphSEIR or running the scripts.

```
julia
]add Plots
]add PyPlot
```

## Installation

To install third-party dependencies and initialize the repository, from the root, simply do:
```
julia
]activate .
]instantiate
```

## Disclaimer

DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

This material is based upon work supported by the Under Secretary of Defense for Research and Engineering under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the Under Secretary of Defense for Research and Engineering.

**© 2022 Massachusetts Institute of Technology.**

Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)

The software/firmware is provided to you on an As-Is basis

Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
