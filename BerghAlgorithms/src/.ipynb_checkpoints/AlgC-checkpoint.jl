using Oscar
using Polymake
using InvertedIndices
using Combinatorics
using LinearAlgebra
export remove!
export getIndex
export isIndependent
export independencyIndex
export isRelevant
export toroidalIndex
export divisorialIndex
export coneContains
export minMaxDivisorial
export BerghC



"""
    remove!(::Array{Int64},::Int64)

    In-place removes a given item from a vector.

# Examples
```jldoctest AlgC
julia> A=[1,2,3,4];

julia> remove!(A,1);

julia> A
[2,3,4]
```
"""
function remove!(a::Array{Int64,1}, item::Int64)
    return deleteat!(a, findall(x->x==item, a))
end

"""
    getIndex(::Array{Int64,1},::Array{Int64,2})

    Returns the first index at which a vector appears as a row of a matrix.

# Examples
```jldoctest AlgC
getIndex([0,1,0],[1 0 0; 0 1 0; 0 0 1])
2
``` 
"""
function getIndex(ray::Array{Int64,1},rayMatrix::Array{Int64,2})
    slice=slicematrix(rayMatrix)
    index=findall(x->x==ray,slice)
    return index[1]
end
    
"""
    isIndependent(::Int64,::Array{Int64,1},::Array{Int64,2})

    Takes a ray matrix, a list of indices representing a cone, and an index represeting a ray of that cone. Determines whether the given ray is independent in the cone (i.e. does not contribute to the multiplicity of the cone).

# Examples
```jldoctest AlgC
julia> isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 2 3])
false

julia> isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 1 1])
true
```
"""
function isIndependent(rayIndex::Int64,cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    if size(cone,1)==1
        return true
    end
    scone=copy(cone)
    subcone=remove!(scone,rayIndex)
    mult=getMultiplicity(cone,rayMatrix)
    submult=getMultiplicity(subcone,rayMatrix)
    if mult==submult
        return true
    else
        return false
    end
end
    
"""
    independencyIndex(::Array{Int64,1},::Array{Int64,2})

    Returns the number of non-independent rays in a cone. Input in indices-ray matrix format.

# Examples
```jldoctest
julia> independencyIndex([1,2,3],[1 0 0 ; 1 2 0; 2 0 3; 0 0 5])
2
```
"""
function independencyIndex(cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    index=0
    for elt in cone
        if isIndependent(elt,cone,rayMatrix)==false
            index+=1
        end
    end
    return index
end
    
"""
    isRelevant(::Array{Int64,1},::Array{Int64,1},::StackyFan)

    Determines if the given ray, relative to the given cone, either has a stacky value greater than 1 or is not independent.

# Examples
```jldoctest AlgC
julia> F=makeStackyFan([1 0 0; 1 2 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([1,2,0],[1,2,3],F)
true

julia> F=makeStackyFan([1 0 0; 0 1 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([0,1,0],[1,2,3],F)
false
"""
function isRelevant(ray::Array{Int64,1},cone::Array{Int64,1},F::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    rayStack=F.stacks[encode(ray)]
    rayIndex=getIndex(ray,rayMatrix)
    rayIndependent=isIndependent(rayIndex,cone,rayMatrix)
    if rayStack != 1 || rayIndependent == false
        return true
    else
        return false
    end
end
    
"""
    toroidalIndex(::Array{Int64,1},::StackyFan,::Dict)
    
    Calculates the toroidal index of the given cone of a divisorial stacky fan, or the number of relevant non-divisorial rays. Compare to divisorialIndex().
    
# Examples
julia> F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);

julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    
julia> toroidalIndex([1,2,3],F,div)
3
julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    
julia> toroidalIndex([1,2,3],F,div)
2
"""
function toroidalIndex(cone::Array{Int64,1},F::StackyFan,div::Dict)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    slice=slicematrix(rayMatrix)
    s=count(x->div[slice[x]]==1,cone)
    flipt=0
    for i in cone
        if div[slice[i]]==0
            if isRelevant(slice[i],cone,F)==false
                flipt+=1
            end
        end
    end
    t=size(cone,1)-flipt
    return t-s
end
    
"""
    divisorialIndex(::Array{Int64,1},::StackyFan,::Dict)

    Calculates the divisorial index (defined by Daniel Bergh) of a given cone in a fan with divisorial rays. Specifically, takes the subcone consisting of all relevant non-divisorial rays in a cone, and counts the number of rays that are relevant in that subcone.

# Examples
julia> F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);

julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    
julia> divisorialIndex([1,2,3],F,div)
3
julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    
julia> divisorialIndex([1,2,3],F,div)
1
"""
function divisorialIndex(cone::Array{Int64,1},F::StackyFan,div::Dict)
    slicedRayMatrix=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS))))
    relRes=Array{Int64,1}[]
    relResStack=Int64[]
    c=0
    for i in cone
        ray=slicedRayMatrix[i]
        stack=F.stacks[encode(ray)]
        if div[ray]==0 && isRelevant(ray,cone,F)==true
            c+=1
            push!(relRes,ray)
            push!(relResStack,stack)
        end
    end
    if c==0
        return 0
    else
        relResIndz=[[i-1 for i in 1:c]]
        relResInd=[i for i in 1:c]
        relResCat=Array{Int64}(transpose(hcat(relRes...)))
        subfan=makeStackyFan(relResCat,relResIndz,relResStack)
        divInd=0
        for ray in relRes
            if isRelevant(ray,relResInd,subfan)==true
                divInd+=1
            end
        end
        return divInd
    end
end
    
"""
    coneContains(::Array{Int64,1},::Array{Int64,1})
    
    Checks whether every index in the first input is also contained in the second input. 
    
# Examples
julia> coneContains([1,2,3],[1,2,3,4])
true
julia> coneContains([1,2,5],[1,2,3,4])
false
"""
function coneContains(A::Array{Int64,1},B::Array{Int64,1})
    out=true
    for i in A
        if !(i in B)
            out=false
        end
    end
    return out
end
    
"""
    minMaxDivisorial(::StackyFan,::Dict)
    
    Calculates the maximal divisorial index of all cones in a stacky fan. Each maximal cone of the fan will contain at most one minimal subcone of maximal divisorial index; a list of such cones is returned.

# Examples
```jldoctest AlgC
julia> F=makeStackyFan([1 2 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);
    
julia> div=Dict([1,2,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    
julia> minMaxDivisorial(F,div)
[[3]]
    
julia> F=makeStackyFan([1 1 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);

julia> div=Dict([1,1,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    
julia> minMaxDivisorial(F,div)
[[1,2,3]]
    
""" 
function minMaxDivisorial(F::StackyFan,div::Dict)
    divMax=0
    coneList=getCones(F.fan)
    divisorialDict=Dict()
    for cone in coneList
        d=divisorialIndex(cone,F,div)
        divisorialDict[cone]=d
        if d>divMax
            divMax=d
        end
    end
    if divMax==0
        return nothing
    end
    divMaxCones=Array{Int64,1}[]
    for cone in coneList
        if divisorialDict[cone]==divMax
            push!(divMaxCones,cone)
        end
    end
    divMaxConesRefined=Array{Int64,1}[]
    maxconeList=convertIncidenceMatrix(F.fan.MAXIMAL_CONES)
    for maxcone in maxconeList
        if divisorialDict[maxcone]==divMax
            maxconeContains=Array{Int64,1}[]
            mincone=maxcone
            for cone in divMaxCones
                if coneContains(cone,maxcone)==true
                    if size(cone,1)<size(mincone,1)
                        mincone=cone
                    end
                end
            end
            if !(mincone in divMaxConesRefined)
                push!(divMaxConesRefined,mincone)
            end
        end
    end
    return divMaxConesRefined
end
    
"""
    BerghC(::StackyFan,::Array{Int64,1})
    
    Takes a stacky fan and a binary array indicating which rays are divisorial, and runs Daniel Bergh's algorithm C. This algorithm performs a series of stacky blowups to reduce the maximal divisorial index of the fan, and returns a fan with a maximal divisorial index of 0.

# Examples
```jldoctest AlgC
julia> F=makeStackyFan([1 1 0;1 3 0; 0 0 1],[[0,1,2]],[1,1,5]);
    
julia> H, div = BerghC(F,[0,0,0]);
    
julia> convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))
5×3 Matrix{Int64}:
 1  1  0
 0  0  1
 2  4  5
 1  3  0
 1  2  0

julia> div
Dict{Any, Any} with 5 entries:
  [1, 1, 0] => 0
  [0, 0, 1] => 1
  [2, 4, 5] => 1
  [1, 2, 0] => 1
  [1, 3, 0] => 0

julia> F=makeStackyFan([1 0;1 3; 5 17],[[0,1],[1,2]],[1,1,5]);
""" 
function BerghC(F::StackyFan,divlist::Array{Int64,1})
    X=deepcopy(F)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
    slicedRayMatrix=slicematrix(rayMatrix)
    div=Dict()
    for i in 1:size(slicedRayMatrix,1)
        div[slicedRayMatrix[i]]=divlist[i]
    end
    while(true)
        rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
        slicedRayMatrix=slicematrix(rayMatrix)
        subdivTargetCones=minMaxDivisorial(X,div)
        if subdivTargetCones==nothing
            break
        end
        blowupList=Array{Array{Int64,1},1}[]
        for cone in subdivTargetCones
            push!(blowupList,slicematrix(rowMinors(rayMatrix,cone)))
        end
        for raycone in blowupList
            indices=Int64[]
            slicedRayMatrix=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS))))
            for ray in raycone
                push!(indices,findall(x->x==ray,slicedRayMatrix)[1])
            end
            cone=sort(indices)
            if size(cone,1)==1
                div[slicedRayMatrix[cone[1]]]=1
            else
                exceptional=findStackyBarycenter(cone,X)
                X=stackyBlowup(X,[x-1 for x in cone], exceptional)
                primExcep=Array{Int64,1}(Polymake.common.primitive(exceptional))
                div[primExcep]=1
            end
        end
    end
    return X, div
end

"""
    BerghC(::StackyFan,::Dict)

    Functions identically to the previous BerghC method, but takes a dictionary that assigns 0 or 1 to the rays of the input fan to indicate divisoriality. In the previous BerghC method, this dictionary is computed.

"""
function BerghC(F::StackyFan,div::Dict)
    X=deepcopy(F)
    while(true)
        rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
        slicedRayMatrix=slicematrix(rayMatrix)
        subdivTargetCones=minMaxDivisorial(F,div)
        if subdivTargetCones==nothing
            break
        end
        blowupList=Set[]
        for cone in subdivTargetCones
            push!(blowupList,slicematrix(rowMinors(rayMatrix,cone)))
        end
        for raycone in blowupList
            indices=Int64[]
            slicedRayMatrix=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS))))
            for ray in raycone
                push!(indices,findall(x->x==ray,slicedRayMatrix)[1])
            end
            cone=sort(indices)
            if size(cone,1)==1
                div[slicedRayMatrix[cone[1]]]=1
            else
                exceptional=findStackyBarycenter(cone,X)
                X=stackyBlowup(X,[x-1 for x in cone], exceptional) 
                div[exceptional]=1
            end
        end
    end
    return X, div
end
