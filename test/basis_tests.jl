@testset "TensorBasis" begin
    for L in 1:4
        basis = TensorBasis(L)
        for i in eachindex(basis)
            @test getposition(basis, getstate(basis, i)) == i
        end
    end

end
