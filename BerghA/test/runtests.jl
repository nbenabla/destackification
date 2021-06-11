using BerghA
using Polymake
using Test


@testset "BerghA.jl" begin
    @testset "test one" begin
        @test 1+1 == 2
    end
    @testset "test two" begin
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);
        F=addStackStructure(X,[1,1]);
        
        @test BerghAEfficient(F,[1,0]).stacks["1,2"] == 5
        
        #New test here
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]])
        F=addStackStructure(X,[1,1])
        @test stackyWeights(BerghAEfficient(F,[1,1])) == [ 5 ,  2 ,  5 ,  10 ]
    end
end
