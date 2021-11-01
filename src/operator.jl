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

abstract type AbstractOperator end
struct SingleBodyOperator{T} <: AbstractOperator 
    site::Int
end

#=
abstract type SigmaX <: AbstractOperator end
abstract type SigmaY <: AbstractOperator end
abstract type SigmaZ <: AbstractOperator end
abstract type SigmaPlus <: AbstractOperator end
abstract type SigmaMinus <: AbstractOperator end
=#

#struct Operator{T} end
#const SigmaX = Operator{:X}()
#const SigmaY = Operator{:Y}()
#const SigmaZ = Operator{:Z}()
#const SigmaPlus = Operator{:+}()
#const SigmaMinus = Operator{:-}()
SigmaX(i::Integer) = SingleBodyOperator{:X}(i)
SigmaY(i::Integer) = SingleBodyOperator{:Y}(i)
SigmaZ(i::Integer) = SingleBodyOperator{:Z}(i)
SigmaPlus(i::Integer) = SingleBodyOperator{:+}(i)
SigmaMinus(i::Integer) = SingleBodyOperator{:-}(i)

eltype(::SingleBodyOperator{:Y}) = ComplexF64




function apply!(state::BitVector, op::SingleBodyOperator{:X}, T::Type=Float64)
    site = op.site
    state[site] = !state[site]
    return one(T)
end
function apply!(state::BitVector, op::SingleBodyOperator{:Y}, T::Type=Float64)
    site = op.site
    s = state[site]
    state[site] = !s
    return s ? -im : im
end
function apply!(state::BitVector, op::SingleBodyOperator{:Z}, T::Type=Float64)
    site = op.site
    s = state[site]
    return s ? -one(T) : one(T)
end

function apply!(state::BitVector, op::SingleBodyOperator{:+}, T::Type=Float64)
    site = op.site
    s = state[site]
    state[site] = !state[site]
    return s ? zero(T) : one(T)
end
function apply!(state::BitVector, op::SingleBodyOperator{:-}, T::Type=Float64)
    site = op.site
    s = state[site]
    state[site] = !state[site]
    return s ? one(T) : zero(T)
end

spmatrix(::SingleBodyOperator{:X}, T::Type=Float64) = sparse([zero(T) one(T); one(T) zero(T)])
spmatrix(::SingleBodyOperator{:Y}, T::Type=ComplexF64) = sparse([zero(T) -im; im zero(T)])
spmatrix(::SingleBodyOperator{:Z}, T::Type=Float64) = sparse([one(T) zero(T); zero(T) one(T)])
spmatrix(::SingleBodyOperator{:+}, T::Type=Float64) = sparse([zero(T) zero(T); one(T) zero(T)])
spmatrix(::SingleBodyOperator{:-}, T::Type=Float64) = sparse([zero(T) one(T); zero(T) zero(T)])

function spmatrix(op::AbstractOperator, basis::AbstractBasis, T::Type=Float64)
    spmatrix(op, basis, basis, T)
end

function spmatrix(op::AbstractOperator, 
    basis1::AbstractBasis, basis2::AbstractBasis,
    T::Type=Float64)

    basis1.L == basis2.L || error("basis1 and basis2 must have same number of sites. Got $(basis1.L) and $(basis2.L)")
    L = basis1.L

    basis_element = BitVector(undef, L)

    rows = Int[]
    cols = Int[]
    values = []

    for i in eachindex(basis1)
        getstate!(basis_element, basis1, i)

        #basis_element[site_index] = !basis_element[site_index]
        val = apply!(basis_element, op, T)

        if basis_element in basis2

            # get the index of the new basis element
            j = getposition(basis2, basis_element)

            # store the indices
            push!(rows, j)
            push!(cols, i)

            # calculate and store the value
            push!(values, val)

        end
    end

    # create the sparse matrices from the rows/columns/values
    spm = sparse(rows, cols, values, length(basis2), length(basis1))

    return spm
end

####################################################################
#
# TensorBasis specialized functions
#
####################################################################

function spmatrix(op::SingleBodyOperator, basis::TensorBasis, T::Type=Float64)
    L = basis.L
    site_index = op.site
    1 <= site_index <= L || error("site must be in [1,$L], got $i")

    m = spmatrix(op, T)
    kron(I(2^(site_index-1)), m, I(2^(L-site_index)))
end


####################################################################
#
# σ hopping models
#
####################################################################

struct SigmaHoppingOBC{T<:Number}
    p::T
    q::T
    α::T
    β::T
    γ::T
    δ::T
end

function SigmaHoppingOBC(args...)
    args = promote(args)
    T = eltype(args)
    SigmaHoppingOBC{T}(args)
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

    L = basis.L

    p = op.p
    q = op.q
    α = op.α
    β = op.β
    γ = op.γ
    δ = op.δ

    D = length(basis)
    H = spzeros(T, D,D)
    H += α * spmatrix(SigmaPlus(1), basis, T)
    H += γ * spmatrix(SigmaMinus(1), basis, T)
    for i in 1:(L-1)
        H += p * spmatrix(SigmaMinus(i), basis, T) * spmatrix(SigmaPlus(i+1), basis, T)
        H += q * spmatrix(SigmaPlus(i), basis, T) * spmatrix(SigmaMinus(i+1), basis, T)
    end
    H += β * spmatrix(SigmaMinus(L), basis, T)
    H += δ * spmatrix(SigmaPlus(L), basis, T)
    H
end