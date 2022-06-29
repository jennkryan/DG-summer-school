push!(LOAD_PATH, ".")
using Pkg
Pkg.add("LinearAlgebra")
Pkg.add(["PyCall","PyPlot"])
Pkg.add(["LinearAlgebra","Einsum", "DelimitedFiles","Delaunay",
"Distributions", "Printf", "CircularArrays","BSplines", "FastGaussQuadrature",
"SpecialPolynomials", "Polynomials", "WriteVTK", "VTKDataIO"])
