module SpinQuantum

using LinearAlgebra
using SparseArrays

export TensorBasis
export getstate, getstate!, getposition

export SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
export spmatrix

export SigmaHoppingOBC

include("basis.jl")
include("operator.jl")







end
