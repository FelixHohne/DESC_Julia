module DESC

using PyCall
using HDF5
using LinearAlgebra
using StaticArrays
using StructArrays
using Interpolations
using HCubature

using PlasmaEquilibriumToolkit

# Write your package code here.

include("descOptimizers.jl")
include("descEquilibrium.jl")
include("descGrid.jl")
include("descGeometry.jl")
include("descObjectives.jl")
include("descContinuation.jl")
include("descLinearObjectives.jl")
include("descUtils.jl")
include("descBasis.jl")
include("descProfiles.jl")
include("descPlot.jl")
include("descPerturbations.jl")
include("descExamples.jl")
include("descIO.jl")

end
