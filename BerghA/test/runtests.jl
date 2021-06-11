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
        
#          @test BerghAEfficient(F,[1,0]) == (StackyFan(Polymake.BigObjectAllocated(Ptr{Nothing} @0x00007f891d77b310), Dict("1,2" => 5, "2,5" => 1, "2,3" => 10, "1,0" => 2, "3,5" => 5, "1,1" => 5)))
        @test BerghAEfficient(F,[1,0]).stacks["1,2"] == 5


    end
end
