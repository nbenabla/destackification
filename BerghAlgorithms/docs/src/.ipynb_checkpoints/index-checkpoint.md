# Destackification

This package implements toric stacky fans in Julia by extending the Oscar/Polymake polyhedral fan functionality. In addition to a variety of helper functions for working with stacky fans, several [algorithms](https://arxiv.org/abs/1409.5713) created by Daniel Bergh are implemented.

## Key Features
- A simple julia struct for polygonal stacky fans, expanding the Oscar/Polymake polyhedral fan struct.
- Helper functions for performing basic operations on stacky fans.
- Algorithms for blowups, simplicialization, and resolution of singularities on fans, adapted from [Macaulay2](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.14/share/doc/Macaulay2/NormalToricVarieties/html/___Normal__Toric__Variety.html).
- Implementation of Bergh's destackification algorithms.
- Visualization of stacky fans in two and three dimensions.
- ???

## Main Functions
```@docs
BerghA
```
