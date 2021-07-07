using BerghAlgorithms
using Polymake
using Test


 @testset "StackyFan Helper Functions" begin
    @test stackyWeights(makeStackyFan([1 0; 1 1; 1 2],[[0,1],[1,2]],[2,2,2]))==[ 2 ,  2 ,  2 ]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 1; 1 2],INPUT_CONES=[[0,1],[1,2]])
    @test stackyWeights(addStackStructure(X,[2,2,2]))==[ 2 ,  2 ,  2 ]

    @test encode([1,0,2,5])=="1,0,2,5"

    F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,2,3,4]);
    @test stackyWeights(F)==[ 1 ,  2 ,  3 ,  4 ]

    F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,2,3,4]);
    @test getRayStack(F,[1,2])==3

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);
    SX = addStackStructure(X, [2,3,5,7]);
    @test stackyWeights(rootConstruction(SX, [1, 4, 2, 1]))==[ 2 ,  12 ,  10 ,  7 ]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);
    SX = StackyFan(X, [2,3,5,7]);
    @test stackyWeights(rootConstructionDistinguishedIndices(SX, [0,1,1,0], [4, 2, 1, 3]))==[ 2 ,  6 ,  5 ,  7 ]
 
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);
    F = StackyFan(X, [2,3,5,7]);
    distinguished = X.RAYS[[2,3],:];
    @test stackyWeights(rootConstructionDistinguished(F, distinguished, [4, 2]))==[ 2 ,  12 ,  10 ,  7 ]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);
    s=[1,2];
    @test convert(Array{Int64,1},vec(Array(findBarycenter(s,X))))==[2,1]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);
    F=addStackStructure(X,[2,3]);
    @test convert(Array{Int64,1},Array(findStackyBarycenter([1,2],F)))==[5,6]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 3; 1 5],INPUT_CONES=[[0,1],[1,2]]);
    F=addStackStructure(X,[1,2,3]);
    @test findStackyRayMatrix(F)==[1 0;4 6;3 15]

    B=[true,true,false,true]
    @test convertBool(B)==[0, 1, 3]

    v=[1,2]
    M=[0 1; 1 1; 1 0]
    @test getConeRank(v,M)==2

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 2 0;0 0 1;0 1 0; 1 1 1],INPUT_CONES=[[0,1,2],[0,2,3,4]])
    @test getDimension(X)==3

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 1],INPUT_CONES=[[0,1,2,3]])
    @test getConeFaces(X,[1,2,3,4],Array(X.RAYS))==[[ 1 ,  2 ], [ 1 ,  3 ], [ 3 ,  4 ], [ 2 ,  4 ]]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])
    @test X.SIMPLICIAL==false
    @test makeSimplicial(X).SIMPLICIAL==true
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0;0 1 0 0;0 0 1 0;1 -1 1 0; 1 0 -2 0],INPUT_CONES=[[0,1,2,3],[0,4]])
    @test X.SIMPLICIAL==false
    @test makeSimplicial(X).SIMPLICIAL==true

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);
    F=addStackStructure(X,[2,3]);
    @test stackyWeights(stackyBlowup(F,[0,1],[1,1]))==[2,1,3]

    F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,1,1,1])
    @test getConesPolymake(F)[1].RAYS==[1 0; 1 1]

    A=[1 2; 3 4];
    @test slicematrix(A)==[[1,2],[3,4]]

    A=[1 2 3;4 5 6; 7 8 9];
    S=Set([1,3]);
    @test rowMinors(A,S)==[1 2 3; 7 8 9]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);
    M=X.MAXIMAL_CONES;
    @test convertIncidenceMatrix(M)==[[1,2],[2,3]]

    C=Polymake.polytope.Cone(INPUT_RAYS=[1 0; 1 2]);
    @test coneMultiplicity(C)==2

    @test getMultiplicity([1,2],[1 0; 1 2; 1 3])==2   

    @test typeof(coneConvert([1, 2, 4],[1 0 0; 0 1 0; 0 0 1; 1 1 1]))==Polymake.BigObjectAllocated

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 1],INPUT_CONES=[[0,1,2,3]]);
    @test Set(getCones(X))==Set([[1,2,3,4],[1,2],[1,3],[3,4],[2,4],[1],[2],[3],[4]])

    C=Polymake.polytope.Cone(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1]);
    @test convert(Array{Int64,2},Array(findFaceContainingRay(C,[1,1,1]).RAYS))==[1 0 1; 1 1 1]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0; 1 0 1 0; 1 1 0 0; 1 1 1 0; 0 0 0 1],INPUT_CONES=[[0,1,2,3,4]]);
    @test findMinimalCone(X,[2,1,0,0])==[1,3]

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0; 1 0 1 0; 1 1 0 0; 1 1 1 0; 0 0 0 1],INPUT_CONES=[[0,1,2,3,4]]);
    @test convert(Array{Int64,2},Array(Polymake.common.primitive(starSubdivision(X,[2,1,0,0]).RAYS)))==[1 0 0 0;1 0 1 0;0 0 0 1;2 1 0 0; 1 1 0 0; 1 1 1 0]

    @test distinguishedAndIntPoint([1,2,4],[1 0 0; 1 2 0;2 1 3; 1 0 3],[1,0,0,0])==true

    @test convertToIncidence([2,3,5],6)==[0,1,1,0,1,0]

    @test compareCones([1,2],[2,3],[1 0 0; 0 1 0; 0 0 1],[1,1,0])==-1
    @test compareCones([1,2],[1,3],[1 0;1 2;1 -1],[1,1,1])==1

    @test extremalCones([[1,2],[2,3],[3,4]],[1 0;1 2; 1 5; 1 8],[0,1,1,0])==[[3,4]]

    C=Polymake.polytope.Cone(INPUT_RAYS=[1 2; 2 1]);
    @test interiorPoints(C)==[[1,1]]

    A=[[1,1,1],[2,1,3],[0,5,4]];
    @test minimalByLex(A)==[0,5,4]

    @test minimalByDist([[0,1,5,7],[3,3,2,2],[8,5,3,6],[2,1,1,10]],[0,1,1,0])==[2,1,1,10]

    @test coneRayDecomposition([1,2,3],[3 5 7; 8 16 9;2 1 3;1 1 1],[2,2,3],[1,1,1,1])==[6,5,52,0]

end

@testset "Visualization" begin
    plot3dSimpCone([[1,0,0],[0,1,0],[0,0,1]]);

    plot2dCone([[1,0],[0,1]]);

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);
    showSimpFan(X);

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);
    F=addStackStructure(X,[1,2]);
    showSimpStackyFan(F);

    C=Polymake.polytope.Cone(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1]);
    @test coneVectorOrder(C)==[[1,1,0],[1,1,1],[1,0,1],[1,0,0]]

    plot3dCone([[1,1,0],[1,1,1],[1,0,1],[1,0,0]]);

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1],INPUT_CONES=[[0,1,2,3]]);
    showFan(X);

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1],INPUT_CONES=[[0,1,2,3]]);
    F=addStackStructure(X,[1,2,2,5]);
    showStackyFan(F);
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
    
# # @testset "slow test" begin
# #          X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[4 1; 7 9],INPUT_CONES=[[0,1]])
# #         F=addStackStructure(X,[1,1])
# #         @test stackyWeights(BerghA(F,[1,0])) == [ 609 ,  29 ,  1 ,  174 ,  29 ,  1740 ,  1218 ,  58 ,  145 ,  1044 ,  1044 ,  290 ,  290 ,  145 ,  406 ,  348 ,  261 ,  14616 ,  14616 ,  609 ,  3480 ,  870 ,  609 ,  174 ,  609 ,  9744 ,  14616 ,  1218 ,  58 ,  145 ,  435 ,  725 ,  1305 ,  1740 ,  6960 ,  3480 ,  870 ,  3480 ,  58464 ,  6090 ,  3480 ,  1392 ,  696 ,  1044 ,  2088 ,  261 ,  174 ,  261 ,  609 ,  406 ,  609 ,  609 ,  406 ,  609 ,  1218 ,  812 ,  1218 ,  2088 ,  261 ,  1044 ,  1160 ,  1740 ,  1740 ,  870 ,  1305 ,  1305 ,  14616 ,  10440 ,  145 ,  435 ,  609 ,  116 ,  580 ,  290 ,  580 ,  1740 ,  3480 ,  3480 ,  261 ,  522 ,  261 ,  522 ,  522 ,  1218 ,  1218 ,  1218 ,  1218 ,  2436 ,  2436 ,  4176 ,  2088 ,  3480 ,  2610 ,  29232 ,  6090 ,  1218 ,  290 ,  145 ,  290 ,  12180 ,  261 ,  522 ]
# # end


@testset "makeSmooth tests" begin
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[4 -1; 0 1],INPUT_CONES=[[0, 1]])
    @test X.SMOOTH_FAN == false
    @test makeSmooth(X).SMOOTH_FAN == true
    @test makeSmooth(X).INPUT_RAYS == [1 -1/4; 0 1; 1 0]
    

    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])
    @test X.SIMPLICIAL == false
    @test X.SMOOTH_FAN == false
    @test makeSmooth(X).SMOOTH_FAN == true
    
    #This test fails, the last ray is given as 4 2 2
    @test makeSmooth(X).INPUT_RAYS == [1 0 0;1 1 0;1 0 1;1 1 1;2 1 1]
    @test makeSmooth(X).SIMPLICIAL == true
    
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0;0 1 0 0;0 0 1 0;1 -1 1 0;1 0 -2 0],INPUT_CONES=[[0,1,2,3],[0,4]])
    @test X.SIMPLICIAL == false
    @test X.SMOOTH_FAN == false
    @test makeSmooth(X).SMOOTH_FAN == true
    @test makeSmooth(X).SIMPLICIAL == true

    
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 2 3;-1 1 1;1 -1 1;-1 -1 1;1 1 -1;-1 1 -1;1 -1 -1;-1 -1 -1],INPUT_CONES=[[0,1,2,3],[0,1,4,5],[0,2,4,6],[1,3,5,7],[2,3,6,7],[4,5,6,7]])
    @test X.SIMPLICIAL == false
    @test X.SMOOTH_FAN == false
    @test makeSmooth(X).SMOOTH_FAN == true
    @test makeSmooth(X).SIMPLICIAL == true
    
    X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;0 1 0;0 0 1;0 -1 -1; -1 0 -1; -2 -1 0],INPUT_CONES=[[0,1,2],[0,1,3],[1,3,4],[1,2,4],[2,4,5],[0,2,5],[0,3,5],[3,4,5]])
    @test X.SIMPLICIAL == true
    @test X.SMOOTH_FAN == false
    @test makeSmooth(X).SMOOTH_FAN == true
    @test makeSmooth(X).SIMPLICIAL == true

end

@testset "toric_blowup tests" begin
       X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])
       @test X.SIMPLICIAL == false
       A = toric_blowup([0,1,2,3], X, nothing)
       @test A.INPUT_RAYS == [1 0 0; 1 1 0; 1 0 1; 1 1 1;2 1 1]
       @test Set(convertIncidenceMatrix(A.INPUT_CONES)) == Set([[1, 2, 5], [1, 3, 5], [3, 4, 5], [2, 4, 5]])
       @test A.SIMPLICIAL == true
       @test A.SMOOTH_FAN == true
    
       A = toric_blowup([0,1], X, nothing)
       @test A.INPUT_RAYS == [1 0 0; 1 1 0; 1 0 1;1 1 1;2 1 0]
       @test Set(convertIncidenceMatrix(A.INPUT_CONES)) == Set([[1, 3, 5], [2, 4, 5], [3, 4, 5]])
       @test A.SIMPLICIAL == true
        
       A = toric_blowup([0,1,2,3], X, [5,3,4])
       @test A.INPUT_RAYS == [1 0 0;1 1 0;1 0 1;1 1 1;5 3 4]
       @test Set(convertIncidenceMatrix(A.INPUT_CONES)) == Set([[1, 2, 5], [1, 3, 5], [2, 4, 5], [3, 4, 5]])
       @test A.SIMPLICIAL == true
    
        A = toric_blowup([0], X, nothing)
        @test A.INPUT_RAYS == [1 0 0;1 1 0;1 0 1;1 1 1]
        @test Set(convertIncidenceMatrix(A.INPUT_CONES)) == Set([[1, 2, 4], [1, 3, 4]])
        @test A.SIMPLICIAL == true
    

       X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]])
       B=toric_blowup([0,1],X,[1,1]);
       @test B.INPUT_RAYS == [1 0;1 2; 1 1]
       @test Set(convertIncidenceMatrix(B.INPUT_CONES)) == Set([[1, 3], [2, 3]])
       @test B.SIMPLICIAL == true
       @test B.SMOOTH_FAN == true
       
end

@testset "AlgC.jl" begin
    A=[1,2,3,4];
    remove!(A,1);
    @test A==[2,3,4]

    @test getIndex([0,1,0],[1 0 0; 0 1 0; 0 0 1])==2

    @test isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 2 3])==false
    @test isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 1 1])==true

    @test independencyIndex([1,2,3],[1 0 0 ; 1 2 0; 2 0 3; 0 0 5])==3

    F=makeStackyFan([1 0 0; 1 2 0; 0 0 1],[[0,1,2]],[1,1,2]);
    @test isRelevant([1,2,0],[1,2,3],F)==true
    F=makeStackyFan([1 0 0; 0 1 0; 0 0 1],[[0,1,2]],[1,1,2]);
    @test isRelevant([0,1,0],[1,2,3],F)==false

    F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);
    div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    @test toroidalIndex([1,2,3],F,div)==3
    div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    @test toroidalIndex([1,2,3],F,div)==2

    F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);
    div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    @test divisorialIndex([1,2,3],F,div)==3
    div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    @test divisorialIndex([1,2,3],F,div)==1
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))== [1 0 0; 0 1 0; 1 1 1; 1 0 2; 1 0 1]
    @test div == Dict([1, 1, 1] => 1, [1, 0, 0] => 0, [1, 0, 1] => 1, [0, 1, 0] => 1, [1, 0, 2] => 0)
    
    @test coneContains([1,2,3],[1,2,3,4])==true
    @test coneContains([1,2,5],[1,2,3,4])==false

    F=makeStackyFan([1 2 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);
    div=Dict([1,2,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    @test minMaxDivisorial(F,div)==[[3]]
    F=makeStackyFan([1 1 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);
    div=Dict([1,1,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    @test minMaxDivisorial(F,div)==[[1,2,3]]

    F=makeStackyFan([1 1 0;1 3 0; 0 0 1],[[0,1,2]],[1,1,5]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[1 1 0; 0 0 1; 2 4 5; 1 3 0; 1 2 0]
    @test div==Dict([1,1,0]=>0,[0,0,1]=>1,[2,4,5]=>1,[1,2,0]=>1,[1,3,0]=>0)

    F=makeStackyFan([1 0;1 3; 5 17],[[0,1],[1,2]],[1,1,5]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[1 0; 2 3; 1 3; 13 44; 5 17]
    @test div==Dict([1,0]=>0,[13,44]=>1,[1,3]=>0,[2,3]=>1,[5,17]=>1)
    
    F=makeStackyFan([1 0;1 3],[[0,1]],[1,1]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[1 0; 2 3; 1 3]
    @test div==Dict([1,0]=>0,[1,3]=>0,[2,3]=>1)
    
    H, div = BerghC(F, [1,1,1]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[1 0; 1 3]
    @test div==Dict([1,0]=>1,[1,3]=>1)
    
    F=makeStackyFan([1 0;1 3],[[0,1]],[3,3]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[1 0; 2 3; 1 3]
    @test div==Dict([1, 0] => 1, [1, 3] => 1, [2, 3] => 1)
    
    
    F=makeStackyFan([0 0 1;1 0 0; 1 3 0],[[0,1,2]],[1,1,1]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[0 0 1; 1 0 0; 2 3 0; 1 3 0]
    @test div==Dict([0, 0, 1] => 0, [2, 3, 0] => 1, [1, 0, 0] => 0, [1, 3, 0] => 0)
    
    F=makeStackyFan([0 0 1;1 0 0; 1 3 0],[[0,1,2]],[1,2,5]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))==[0 0 1; 1 0 0; 7 15 0; 1 3 0]
    @test div==Dict([7, 15, 0] => 1, [0, 0, 1] => 0, [1, 0, 0] => 1, [1, 3, 0] => 1)
    
    F = makeStackyFan([1 0; 1 1; 1 2],[[0,1],[1,2]],[2,2,2]);
    H, div = BerghC(F,[0,0,0]);
    @test convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS)) ==[1 0; 2 1; 1 1; 2 3; 1 2]
    @test div==Dict([1, 0] => 1, [1, 2] => 1, [2, 1] => 1, [1, 1] => 1, [2, 3] => 1)
    
    
end

