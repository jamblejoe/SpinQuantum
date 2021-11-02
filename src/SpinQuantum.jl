module SpinQuantum

using LinearAlgebra
using SparseArrays

export TensorBasis, TotalSpinConservedBasis
export getstate, getstate!, getposition

export SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
export spmatrix

export SigmaHoppingChainOBC, SigmaHoppingChainPBC

include("basis.jl")
include("operator.jl")


end
