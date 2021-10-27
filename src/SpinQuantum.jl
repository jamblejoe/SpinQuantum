module SpinQuantum

using LinearAlgebra
using SparseArrays

export TensorBasis
export getstate, getstate!, getposition

export SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
export spmatrix
#export σ_x_spmatrix


abstract type AbstractBasis end

struct TensorBasis
    L::Int
end

function getposition(basis::TensorBasis, state::AbstractVector)
    L = basis.L

	state in basis || error("state $state not in basis")
	index = 1
	for i in eachindex(state)
		index += 2^(L-i) * state[i]
	end
	index
end

function getstate(basis::TensorBasis, index::Integer)
    getstate!(BitVector(undef, basis.L), basis, index)
end

function getstate!(state::AbstractVector, basis::TensorBasis, index::Integer)
    L = basis.L

    1<= index <= length(basis) || error("Index $index out of bounds [1,$(length(basis))]")

    digits!(state, index-1, base=2)
    reverse!(state)
    state
end

#sites(basis::TensorBasis) = basis.L

Base.size(basis::TensorBasis) = (length(basis),)
Base.length(basis::TensorBasis) = 2^basis.L
Base.isequal(b1::TensorBasis, b2::TensorBasis) = b1.L == b2.L

# Iterator interface
function Base.iterate(basis::TensorBasis, index::Integer=1)
    1 <= index <= length(basis) ? (getstate(basis, index), index+1) : nothing
end

Base.eltype(::TensorBasis) = BitVector

# Array interface
"""
Returns the basis element at position index.
"""
Base.getindex(basis::TensorBasis, index::Integer) = getstate(basis, index)
Base.firstindex(::TensorBasis) = 1
Base.lastindex(basis::TensorBasis) = length(basis)

Base.in(state::AbstractVector, basis::TensorBasis) = length(state)==basis.L && all(0 .<= state .<= 1)
Base.in(state::BitVector, basis::TensorBasis) = length(state)==basis.L
Base.eachindex(basis::TensorBasis) = 1:length(basis)


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

#abstract type AbstractOperator end
#struct Operator{T} <: AbstractOperator end
#=
abstract type SigmaX <: AbstractOperator end
abstract type SigmaY <: AbstractOperator end
abstract type SigmaZ <: AbstractOperator end
abstract type SigmaPlus <: AbstractOperator end
abstract type SigmaMinus <: AbstractOperator end
=#

struct Operator{T} end
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




end
