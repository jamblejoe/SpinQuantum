module SpinQuantum

using LinearAlgebra
using SparseArrays

export TensorBasis, TotalSpinConservedBasis
export getstate, getstate!, getposition

export SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
export SigmaPlusMinus, SigmaMinusPlus
export spmatrix

include("basis.jl")
include("operator.jl")


end
