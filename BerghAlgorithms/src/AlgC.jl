export isIndependent, isRelevant
export independencyIndex, toroidalIndex, divisorialIndex, minMaxDivisorial
export remove!, getIndex, coneContains
export BerghC


"""
    remove!(::Array{Int64},::Int64)

    In-place removes a given item from a vector.

# Examples
```jldoctest
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
```jldoctest
getIndex([0,1,0],[1 0 0; 0 1 0; 0 0 1])
2
``` 
"""
function getIndex(ray::Array{Int64,1},rayMatrix::Array{Int64,2})
    slice=slicematrix(rayMatrix)
    index=findfirst(x->x==ray,slice)
    return index
end
    
"""
    isIndependent(::Int64,::Array{Int64,1},::Array{Int64,2})

    Takes a ray matrix, a list of indices representing a cone, and an index represeting a ray of that cone. Determines whether the given ray is independent in the cone (i.e. does not contribute to the multiplicity of the cone).

# Examples
```jldoctest
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
    return mult==submult
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
    return count(elt -> !isIndependent(elt, cone, rayMatrix), cone)
end
    
"""
    isRelevant(::Array{Int64,1},::Array{Int64,1},::StackyFan)

    Determines if the given ray, relative to the given cone, either has a stacky value greater than 1 or is not independent.

# Examples
```jldoctest
julia> F=makeStackyFan([1 0 0; 1 2 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([1,2,0],[1,2,3],F)
true

julia> F=makeStackyFan([1 0 0; 0 1 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([0,1,0],[1,2,3],F)
false
```
"""
function isRelevant(ray::Array{Int64,1},cone::Array{Int64,1},F::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    rayStack=F.stacks[encode(ray)]
    rayIndex=getIndex(ray,rayMatrix)
    rayIndependent=isIndependent(rayIndex,cone,rayMatrix)
    return rayStack != 1 || !rayIndependent
end
    
"""
    toroidalIndex(::Array{Int64,1},::StackyFan,::Dict)
    
    Calculates the toroidal index of the given cone of a divisorial stacky fan, or the number of relevant non-divisorial rays. Compare to divisorialIndex().
    
# Examples
```jldoctest
julia> F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);

julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    
julia> toroidalIndex([1,2,3],F,div)
3
julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    
julia> toroidalIndex([1,2,3],F,div)
2
```
"""
function toroidalIndex(cone::Array{Int64,1},F::StackyFan,div::Dict)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    slice=slicematrix(rayMatrix)
    # Find number of divisorial rays 
    # div is a dictionary that represents which rays are divisorial, 0 represents a non-divisorial ray and 1 represents a divisorial ray
    s=count(x->div[slice[x]]==1,cone)
    #flipt counts the number of non-divisorial and irrelevant rays
    flipt = count(i -> div[slice[i]]==0 && !isRelevant(slice[i], cone, F), cone)
    # number of relevant cones
    t=size(cone,1)-flipt
    # return the toroidal index, which is the number of relevant residual (non-divisorial) rays
    return t-s
end
    
"""
    divisorialIndex(::Array{Int64,1},::StackyFan,::Dict)

    Calculates the divisorial index (defined by Daniel Bergh) of a given cone in a fan with divisorial rays. 
    Specifically, takes the subcone consisting of all relevant non-divisorial rays in a cone, and counts the number of rays that are relevant in that subcone.

# Examples
```jldoctest
julia> F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1]);

julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>0);
    
julia> divisorialIndex([1,2,3],F,div)
3
julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1);
    
julia> divisorialIndex([1,2,3],F,div)
1
```
"""
function divisorialIndex(cone::Array{Int64,1},F::StackyFan,div::Dict)
    slicedRayMatrix=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS))))
    # relevant residual cones
    relRes=Array{Int64,1}[]
    relResStack=Int64[]
    # c is the number of non-divisorial relevant cones
    c=0
    for i in cone
        ray=slicedRayMatrix[i]
        stack=F.stacks[encode(ray)]
        # If the ray is non-divisorial and relevant, increment c by one, add ray to relRes, and stack to relResStack. The rays in relRes are used to build a new cone.
        if div[ray]==0 && isRelevant(ray,cone,F)
            c+=1
            push!(relRes,ray)
            push!(relResStack,stack)
        end
    end
    # If there are no relevant residual cones in F, the divisorial index is 0
    if c==0
        return 0
    else
        # convert to 0-indexing
        relResIndz=[[i-1 for i in 1:c]]
        relResInd=[i for i in 1:c]
        relResCat=Array{Int64}(transpose(hcat(relRes...)))
        # Construct a subfan consisting of all relevant residual rays
        subfan=makeStackyFan(relResCat,relResIndz,relResStack)
        divInd=0
        for ray in relRes
            # Iterate over relevant residual cones to count the number of relevant rays in the subfan
            # This count represents the divisorial index
            if isRelevant(ray,relResInd,subfan)
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
```jldoctest
julia> coneContains([1,2,3],[1,2,3,4])
true
julia> coneContains([1,2,5],[1,2,3,4])
false
```
"""
function coneContains(A::Array{Int64,1},B::Array{Int64,1})
    return issubset(A, B)
end
    
"""
    minMaxDivisorial(::StackyFan,::Dict)
    
    Calculates the maximal divisorial index of all cones in a stacky fan. 
    Each maximal cone of the fan will contain at most one minimal subcone of maximal divisorial index; a list of such cones is returned.

# Examples
```jldoctest
julia> F=makeStackyFan([1 2 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);
    
julia> div=Dict([1,2,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    
julia> minMaxDivisorial(F,div)
[[3]]
    
julia> F=makeStackyFan([1 1 0;1 3 0; 3 0 1],[[0,1,2]],[1,1,5]);

julia> div=Dict([1,1,0]=>0,[1,3,0]=>0,[3,0,1]=>0);
    
julia> minMaxDivisorial(F,div)
[[1,2,3]]
```
""" 
function minMaxDivisorial(F::StackyFan,div::Dict)
    # Calculates the maximal divisorial index of any cone in the fan
    divMax=0
    coneList=getCones(F.fan)
    # dictionary that represents each cone with its divisorial index
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

    # cones with maximal divisorial index
    divMaxCones=Array{Int64,1}[]
    for cone in coneList
        if divisorialDict[cone]==divMax
             # if the cone's divisorial index is the fan's maximal divisorial index, add the cone to divMaxCones
            push!(divMaxCones,cone)
        end
    end

    #divMaxConesRefined stores the cones in divMaxCones that are minimal with respect to inclusion
    divMaxConesRefined=Array{Int64,1}[]
    # List of maximal cones in F
    maxconeList=convertIncidenceMatrix(F.fan.MAXIMAL_CONES)
    for maxcone in maxconeList
        # if the div index of the current maxcone is the fan's max div index, its minimal subcone with maximal divisorial index is calculated
        if divisorialDict[maxcone]==divMax
            maxconeContains=Array{Int64,1}[]
            mincone=maxcone
            for cone in divMaxCones
                if coneContains(cone,maxcone) && size(cone,1)<size(mincone,1)
                    mincone=cone
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
    
    Takes a stacky fan and a binary array indicating which rays are divisorial, and runs Daniel Bergh's algorithm C. 
    This algorithm performs a series of stacky blowups to reduce the maximal divisorial index of the fan, and returns a fan with a maximal divisorial index of 0.

# Examples
```jldoctest
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

julia> H, div = BerghC(F,[0,0,0]);

julia> convert(Array{Int64,2},Polymake.common.primitive(H.fan.RAYS))
5×2 Matrix{Int64}:
  1   0
  2   3
  1   3
 13  44
  5  17
julia> div
Dict{Any, Any} with 5 entries:
  [1, 0]   => 0
  [13, 44] => 1
  [1, 3]   => 0
  [2, 3]   => 1
  [5, 17]  => 1
```
""" 
function BerghC(F::StackyFan,divlist::Array{Int64,1})
    X=deepcopy(F)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
    slicedRayMatrix=slicematrix(rayMatrix)
    div=Dict()
    # Populate dictionary of divisorial rays (div) from divlist
    for i in 1:size(slicedRayMatrix,1)
        div[slicedRayMatrix[i]]=divlist[i]
    end
    while(true)
        rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
        slicedRayMatrix=slicematrix(rayMatrix)
        # Find the cones with maximal divisorial index in X
        subdivTargetCones=minMaxDivisorial(X,div)
        # If there are no such cones, the algorithm terminates
        if subdivTargetCones==nothing
            break
        end
        blowupList=Array{Array{Int64,1},1}[]
        # Iterate through cones with maximal divisorial index
        for cone in subdivTargetCones
            # Add each cone's ray representation to blowupList
            push!(blowupList,slicematrix(rowMinors(rayMatrix,cone)))
        end
        for raycone in blowupList
            indices=Int64[]
            slicedRayMatrix=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS))))
            for ray in raycone
                push!(indices,findall(x->x==ray,slicedRayMatrix)[1])
            end
            # Sort indices in ascending order
            cone=sort(indices)
            if size(cone,1)==1
                div[slicedRayMatrix[cone[1]]]=1
            else
                exceptional=findStackyBarycenter(cone,X)
                # perform the blowup
                X=stackyBlowup(X,[x-1 for x in cone], exceptional)
                # convert exceptional ray to its primitive form
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
