## Functions

```@meta StackyFan
DocTestSetup = quote
    using BerghAlgorithms
end
```

## Main Functions for Fans
```@docs
BerghAlgorithms.toric_blowup
BerghAlgorithms.makeSimplicial
BerghAlgorithms.makeSmooth
BerghAlgorithms.starSubdivision
```

## Main Functions for Destackification
```@docs
BerghAlgorithms.BerghA
BerghAlgorithms.stackyBlowup
BerghAlgorithms.BerghC
```

## Main Visualization Functions
```@docs
BerghAlgorithms.showFan
BerghAlgorithms.showStackyFan
```

## Helper Functions for Fans
```@docs
BerghAlgorithms.findBarycenter
BerghAlgorithms.convertBool
BerghAlgorithms.getConeRank
BerghAlgorithms.getDimension
BerghAlgorithms.getConeFaces
BerghAlgorithms.slicematrix
BerghAlgorithms.rowMinors
BerghAlgorithms.convertIncidenceMatrix
BerghAlgorithms.coneMultiplicity
BerghAlgorithms.getMultiplicity
BerghAlgorithms.coneConvert
BerghAlgorithms.getCones
BerghAlgorithms.findFaceContainingRay
BerghAlgorithms.findMinimalCone
BerghAlgorithms.convertToIncidence
BerghAlgorithms.interiorPoints
```

## Helper Functions for Stacky Fans
```@docs
BerghAlgorithms.StackyFan
BerghAlgorithms.makeStackyFan
BerghAlgorithms.addStackStructure
BerghAlgorithms.encode
BerghAlgorithms.stackyWeights
BerghAlgorithms.getRayStack
BerghAlgorithms.rootConstruction
BerghAlgorithms.rootConstructionDistinguished
BerghAlgorithms.rootConstructionDistinguishedIndices
BerghAlgorithms.findStackyBarycenter
BerghAlgorithms.findStackyRayMatrix
BerghAlgorithms.getConesPolymake
BerghAlgorithms.distinguishedAndIntPoint
BerghAlgorithms.compareCones
BerghAlgorithms.extremalCones
BerghAlgorithms.minimalByLex
BerghAlgorithms.minimalByDist
BerghAlgorithms.coneRayDecomposition
BerghAlgorithms.remove!
BerghAlgorithms.getIndex
BerghAlgorithms.isIndependent
BerghAlgorithms.independencyIndex
BerghAlgorithms.isRelevant
BerghAlgorithms.toroidalIndex
BerghAlgorithms.divisorialIndex
BerghAlgorithms.coneContains
BerghAlgorithms.minMaxDivisorial
```

## Visualization Helper Functions
```@docs
BerghAlgorithms.plot3dSimpCone
BerghAlgorithms.plot2dCone
BerghAlgorithms.showSimpFan
BerghAlgorithms.showSimpStackyFan
BerghAlgorithms.coneVectorOrder
BerghAlgorithms.plot3dCone
```