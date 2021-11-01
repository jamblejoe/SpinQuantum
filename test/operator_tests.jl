@testset "Sigma Hopping OBC" begin

    @testset "L = 2" begin
        
        basis = TensorBasis(2)

        for (p,q,α,β,γ,δ) in Iterators.product([[0.14, 0.33, 0.79] for _ in 1:6]...)
            H = spmatrix(SigmaHoppingOBC(p,q,α,β,γ,δ), basis)
            H_correct = [
                0 β γ 0;
                δ 0 p γ;
                α q 0 β;
                0 α δ 0
            ]
            @test H == H_correct
        end
    end

    @testset "L = 3" begin
        
        basis = TensorBasis(3)

        for (p,q,α,β,γ,δ) in Iterators.product([[0.14, 0.33, 0.79] for _ in 1:6]...)
            H = spmatrix(SigmaHoppingOBC(p,q,α,β,γ,δ), basis)
            H_correct = [
            #   1 2 3 4 5 6 7 8
                0 β 0 0 γ 0 0 0; #1
                δ 0 p 0 0 γ 0 0; #2
                0 q 0 β p 0 γ 0; #3
                0 0 δ 0 0 p 0 γ; #4
                α 0 q 0 0 β 0 0; #5
                0 α 0 q δ 0 p 0; #6
                0 0 α 0 0 q 0 β; #7
                0 0 0 α 0 0 δ 0  #8
            ]
            @test H == H_correct
        end
    end

end
