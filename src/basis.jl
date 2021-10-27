abstract type AbstractBasis end

struct TensorBasis <: AbstractBasis
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