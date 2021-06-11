using BerghAlgorithms
using Polymake
using Test


@testset "StackyFan.jl" begin
    @test slicematrix([1 2; 3 4]) == [[ 1 ,  2 ], [ 3 ,  4 ]]
end

@testset "BerghA.jl" begin
    @testset "test stackyWeights" begin
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);
        F=addStackStructure(X,[1,1]);
        @test BerghA(F,[1,0]).stacks["1,2"] == 5
        @test stackyWeights(BerghA(F,[1,1])) == [ 5 ,  2 ,  5 ,  10 ]
        @test stackyWeights(BerghA(F,[1,0])) == [ 5 ,  5 ,  1 ,  2 ,  5 ,  10 ]
        @test stackyWeights(BerghA(F,[0,1])) == [ 1 ,  5 ,  5 ,  2 ,  5 ,  1 ,  5 ,  10 ]
        
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 3; 4 5 6; 2 3 1],INPUT_CONES=[[0,1,2]])
        F=addStackStructure(X,[1,1,1])
        @test stackyWeights(BerghA(F,[1,1,1])) == [ 28 ,  21 ,  84 ,  28 ,  84 ,  84 ,  42 ,  84 ,  168 ]

        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 2; 2 1 1; 5 3 9],INPUT_CONES=[[0,1,2]])
        F=addStackStructure(X,[1,1,1])
        @test stackyWeights(BerghA(F,[1,1,1])) == [ 4 ,  4 ,  8 ,  4 ,  4 ,  8 ]
        
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 4 3; 1 5],INPUT_CONES=[[0,1],[1,2]])
        F=addStackStructure(X,[1,1,1])
        @test stackyWeights(BerghA(F,[1,1,1])) == [ 4 ,  34 ,  6 ,  68 ,  34 ,  68 ,  34 ,  6 ,  12 ]
        
        X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[323 179; 44 135],INPUT_CONES=[[0,1]])
        F=addStackStructure(X,[1,1])
        @test stackyWeights(BerghA(F,[1,1])) == [ 491602456800 ,  49160245680 ,  294961474080 ,  468192816 ,  12173013216 ,  73038079296 ,  12173013216 ,  97384105728 ,  73038079296 ,  196640982720 ,  131093988480 ,  156064272 ,  936385632 ,  468192816 ,  196640982720 ,  49160245680 ,  9832049136 ,  292152317184 ,  5899229481600 ,  3932819654400 ,  589922948160 ,  393281965440 ,  12173013216 ,  8115342144 ,  12173013216 ,  5899229481600 ,  589922948160 ,  737403685200 ,  51126655507200 ,  292152317184 ,  51126655507200 ,  5899229481600 ,  5899229481600 ,  589922948160 ,  589922948160 ,  196640982720 ,  2949614740800 ,  5899229481600 ,  294961474080 ,  589922948160 ,  196640982720 ,  589922948160 ,  393281965440 ,  936385632 ,  24346026432 ,  24346026432 ,  11798458963200 ,  1179845896320 ,  737403685200 ,  1474807370400 ,  2949614740800 ,  5899229481600 ,  589922948160 ,  11798458963200 ,  1179845896320 ,  2949614740800 ,  1474807370400 ,  5899229481600 ,  102253311014400 ]
        
    end
end
    
@testset "slow test" begin
         X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[4 1; 7 9],INPUT_CONES=[[0,1]])
        F=addStackStructure(X,[1,1])
        @test stackyWeights(BerghA(F,[1,0])) == [ 609 ,  29 ,  1 ,  174 ,  29 ,  1740 ,  1218 ,  58 ,  145 ,  1044 ,  1044 ,  290 ,  290 ,  145 ,  406 ,  348 ,  261 ,  14616 ,  14616 ,  609 ,  3480 ,  870 ,  609 ,  174 ,  609 ,  9744 ,  14616 ,  1218 ,  58 ,  145 ,  435 ,  725 ,  1305 ,  1740 ,  6960 ,  3480 ,  870 ,  3480 ,  58464 ,  6090 ,  3480 ,  1392 ,  696 ,  1044 ,  2088 ,  261 ,  174 ,  261 ,  609 ,  406 ,  609 ,  609 ,  406 ,  609 ,  1218 ,  812 ,  1218 ,  2088 ,  261 ,  1044 ,  1160 ,  1740 ,  1740 ,  870 ,  1305 ,  1305 ,  14616 ,  10440 ,  145 ,  435 ,  609 ,  116 ,  580 ,  290 ,  580 ,  1740 ,  3480 ,  3480 ,  261 ,  522 ,  261 ,  522 ,  522 ,  1218 ,  1218 ,  1218 ,  1218 ,  2436 ,  2436 ,  4176 ,  2088 ,  3480 ,  2610 ,  29232 ,  6090 ,  1218 ,  290 ,  145 ,  290 ,  12180 ,  261 ,  522 ]
end
