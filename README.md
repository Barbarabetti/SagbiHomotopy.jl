# SagbiHomotopy

[![Build Status](https://github.com/Barbarabetti/SagbiHomotopy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Barbarabetti/SagbiHomotopy.jl/actions/workflows/CI.yml?query=branch%3Amain)

SagbiHomotopy is a Julia package to solve equations that are expressed as linear combination of a SAGBI basis. A detailed documentation is available at the MathRepo page [MathRepo page SagbiHomotopy](https://mathrepo.mis.mpg.de/SagbiHomotopy/).

## Installation

```julia
using Pkg
Pkg.add("SagbiHomotopy")
```
SagbiHomotopy is supported on Julia 1.11.3 and later.

## Example

The following is a minimum working example.

```julia
using SagbiHomotopy
using HomotopyContinuation
using Oscar 

@var x,y, z

sagbi = [[x,y,(x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
w = get_weight(sagbi)

degree_map(sagbi)
degree_monomial_map(sagbi,w)

@var p[1:4]
@var q[1:4]

lin_sys = [ rand(Int64,2,4)*p , rand(Int64,1,4)*q]

@time sagbi_homotopy(lin_sys, sagbi; weight = w)
```

<!-- ## Citation

If you found this package to be useful in academic work, then please cite:

```bibtex
@article{
}
``` -->