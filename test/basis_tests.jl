@testset "TensorBasis" begin

    @testset "L=2" begin
        
        basis = TensorBasis(2)
        @test getstate(basis, 1) == [0,0]
        @test getstate(basis, 2) == [0,1]
        @test getstate(basis, 3) == [1,0]
        @test getstate(basis, 4) == [1,1]
    end

    @testset "L=3" begin
        
        basis = TensorBasis(3)
        @test getstate(basis, 1) == [0,0,0]
        @test getstate(basis, 2) == [0,0,1]
        @test getstate(basis, 3) == [0,1,0]
        @test getstate(basis, 4) == [0,1,1]
        @test getstate(basis, 5) == [1,0,0]
        @test getstate(basis, 6) == [1,0,1]
        @test getstate(basis, 7) == [1,1,0]
        @test getstate(basis, 8) == [1,1,1]
    end

    for L in 1:6
        basis = TensorBasis(L)
        for i in eachindex(basis)
            @test getposition(basis, getstate(basis, i)) == i
        end
    end

end
