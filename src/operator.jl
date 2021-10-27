"""
Pauli matrices
σ_x = [0 1; 1 0]
σ_y = [0 -im; im 0]
σ_z = [1 0; 0 -1]

σ_+ = [0 0; 1 0] = 1/2 (σ_x - iσ_y)
σ_- = [0 1; 0 0] = 1/2 (σ_x + iσ_y)

|0> = [1, 0]'
|1> = [0, 1]'
"""

function σ_x_spmatrix(site_index::Integer, basis::AbstractBasis, T::Type=Float64)

    L = basis.L
    basis_length = length(basis)

    1 <= site_index <= L || error("site must be in [1,$L], got $i")

    basis_element = BitVector(undef, L)

    rows = Int[]
    cols = Int[]
    values = T[]

    for i in eachindex(basis)
        getstate!(basis_element, basis, i)

        # flip the bit
        basis_element[site_index] = !basis_element[site_index]

        if basis_element in basis

            # get the index of the new basis element
            j = getposition(basis, basis_element)

            # store the indices
            push!(rows, j)
            push!(cols, i)

            # calculate and store the value
            value = one(T)
            push!(values, value)

        end
    end

    # create the sparse matrices from the rows/columns/values
    spm = sparse(rows, cols, values, basis_length, basis_length)

    return spm
end

####################################################################
#
# TensorBasis specialized functions
#
####################################################################

abstract type AbstractOperator end
struct Operator{T} <: AbstractOperator end
#=
abstract type SigmaX <: AbstractOperator end
abstract type SigmaY <: AbstractOperator end
abstract type SigmaZ <: AbstractOperator end
abstract type SigmaPlus <: AbstractOperator end
abstract type SigmaMinus <: AbstractOperator end
=#

#struct Operator{T} end
const SigmaX = Operator{:X}()
const SigmaY = Operator{:Y}()
const SigmaZ = Operator{:Z}()
const SigmaPlus = Operator{:+}()
const SigmaMinus = Operator{:-}()

spmatrix(::Operator{:X}, T::Type=Float64) = sparse([zero(T) one(T); one(T) zero(T)])
spmatrix(::Operator{:Y}, T::Type=ComplexF64) = sparse([zero(T) -im; im zero(T)])
spmatrix(::Operator{:Z}, T::Type=Float64) = sparse([one(T) zero(T); zero(T) one(T)])
spmatrix(::Operator{:+}, T::Type=Float64) = sparse([zero(T) zero(T); one(T) zero(T)])
spmatrix(::Operator{:-}, T::Type=Float64) = sparse([zero(T) one(T); zero(T) zero(T)])



function spmatrix(op::Operator, site_index::Integer, basis::TensorBasis, T::Type=Float64)
    L = basis.L

    1 <= site_index <= L || error("site must be in [1,$L], got $i")

    m = spmatrix(op, T)
    kron(I(2^(site_index-1)), m, I(2^(L-site_index)))
end


####################################################################
#
# σ hopping models
#
####################################################################

struct SigmaHoppingOBC{T<:AbstractFloat} <: AbstractOperator
    p::T
    q::T
    α::T
    β::T
    γ::T
    δ::T
end

function SigmaHoppingOBC(parameters::Dict)
    p = parameters["p"]
    q = parameters["q"]

    α = haskey(parameters, "α") ? parameters["α"] : parameters["alpha"]
    β = haskey(parameters, "β") ? parameters["β"] : parameters["beta"]
    γ = haskey(parameters, "γ") ? parameters["γ"] : parameters["gamma"]
    δ = haskey(parameters, "δ") ? parameters["δ"] : parameters["delta"]

    SigmaHoppingOBC(p,q,α,β,γ,δ)
end

function spmatrix(op::SigmaHoppingOBC,
    basis::AbstractBasis, T::Type=Float64)

    p = op.p
    q = op.q
    α = op.α
    β = op.β
    γ = op.γ
    δ = op.δ

    D = length(basis)
    H = spzeros(T, D,D)
    H += α * spmatrix(SigmaPlus, 1, basis)
    H += γ * spmatrix(SigmaMinus, 1, basis)
    for i in 1:(L-1)
        H += p * spmatrix(SigmaMinus, i, basis) * spmatrix(SigmaPlus, i+1, basis)
        H += q * spmatrix(SigmaPlus, i, basis) * spmatrix(SigmaMinus, i+1, basis)
    end
    H += β * spmatrix(SigmaMinus, L, basis)
    H += δ * spmatrix(SigmaPlus, L, basis)
    H
end