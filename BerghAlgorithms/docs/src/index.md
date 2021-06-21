# Destackification

This package implements toric stacky fans in Julia by extending the Oscar/Polymake polyhedral fan functionality. In addition to a variety of helper functions for working with stacky fans, several [algorithms](https://arxiv.org/abs/1409.5713) created by Daniel Bergh are implemented.

## Key Features
- A simple julia struct for polygonal stacky fans, expanding the Oscar/Polymake polyhedral fan struct.
- Helper functions for performing basic operations on stacky fans.
- Algorithms for blowups, simplicialization, and resolution of singularities on fans, adapted from [Macaulay2](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.14/share/doc/Macaulay2/NormalToricVarieties/html/___Normal__Toric__Variety.html).
- Implementation of Bergh's destackification algorithms.
- Visualization of stacky fans in two and three dimensions.
- ???

## Main Functions for Fans
```@docs
BerghAlgorithms.toric_blowup
BerghAlgorithms.makeSimplicial
BerghAlgorithms.makeSmooth
BerghAlgorithms.starSubdivision
```

## Main Functions for Stacky Fans
```@docs
BerghA
BerghAlgorithms.stackyBlowup
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
BerghAlgorithms.coneConvert
BerghAlgorithms.getCones
BerghAlgorithms.findFaceContainingRay
BerghAlgorithms.findMinimalCone
BerghAlgorithms.convertToIncidence
BerghAlgorithms.interiorPoints
```

## Helper Functions for Stacky Fans
```@docs
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
```

## Visualization Helper Functions
```@docs
BerghAlgorithms.plot3dSimpCone
BerghAlgorithms.plot2dCone
BerghAlgorithms.showSimpFan
BerghAlgorithms.showSimpStackyFan
BerghAlgorithms.coneVectorOrder
BerghAlgorithms.plot3dCone
BerghAlgorithms.showFan
BerghAlgorithms.showStackyFan
```
