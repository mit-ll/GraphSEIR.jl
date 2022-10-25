# © 2022 Massachusetts Institute of Technology.  See LICENSE file for details.

push!(LOAD_PATH, "../src/")
using Documenter, Pkg
Pkg.activate("../")
using GraphSEIR
makedocs(sitename="GraphSEIR Documentation",
        modules=[GraphSEIR]
        #pages=[
        #       "Home" => "index.md"
        #]
)
