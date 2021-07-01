using Oscar
using Polymake
using InvertedIndices
using Combinatorics
using LinearAlgebra
# StackyFan struct
export StackyFan, makeStackyFan, addStackStructure
export stackyWeights, getRayStack, encode
# Fan
export toric_blowup, makeSimplicial, makeSmooth, starSubdivision
export getConeRank, getDimension, getConeFaces
export coneMultiplicity, getMultiplicity
export findFaceContainingRay, findMinimalCone, interiorPoints, findBarycenter
export rootConstruction, rootConstructionDistinguished, rootConstructionDistinguishedIndices
# Stacky fan
export stackyBlowup, coneRayDecomposition
export findStackyBarycenter, findStackyRayMatrix
# Utility
export slicematrix, rowMinors, convertIncidenceMatrix, convertToIncidence
export coneConvert, getCones, getConesPolymake
# Display
export plot3dSimpCone, plot2dCone, showSimpFan, showSimpStackyFan, coneVectorOrder, plot3dCone
export showFan, showStackyFan


"""
    Structure to store information of a stacky fan - this is a fan together with a dictionary assigning stacky values to each ray.

# Properties:
-> `fan` - the underlying fan, as a polymake object
-> `scalars` - the array of stacky values
-> `stacks` - a dictionary assigning stacky values to each ray.

"""
struct StackyFan
    fan::Polymake.BigObjectAllocated
    stacks::Dict{String, Int64}
    # Constructors for the StackyFan object
    StackyFan(
        fan::Polymake.BigObjectAllocated,
        stacks::Dict{String, Int64}) = new(fan, stacks)
    StackyFan(
        rays::Array{Int64, 2},
        cones::Array{Array{Int64, 1}, 1},
        scalars::Array{Int64, 1}) = makeStackyFan(rays, cones, scalars)
    StackyFan(
        fan::Polymake.BigObjectAllocated,
        scalars::Array{Int64, 1}) = addStackStructure(fan, scalars)
end

## Helper functions

"""
    makeStackyFan(::Array{Int64,2},::Array{Array{Int64,1},1},::Array{Int64,1}))

Function to generate a stacky fan from a matrix representing rays as row vectors, a vector of vectors representing the rays contained in each cone, and a vector of stacky values to be assigned the rays. The second input should be zero-indexed.

# Examples
```jldoctest
julia> makeStackyFan([1 0; 1 1; 1 2],[[0,1],[1,2]],[2,2,2])
[ 2 ,  2 ,  2 ]
```
"""
function makeStackyFan(
    rays::Array{<:Number,2},
    cones::Array{Array{Int64,1},1},
    scalars::Array{Int64,1})

    # Construct a normal fan from the given rays and cones
    fan = fulton.NormalToricVariety(INPUT_RAYS=rays, INPUT_CONES=cones)
    
    # Construct the dictionary
    stack_rays = mapslices(encode, convert(Array{Int64,2},Array(Polymake.common.primitive(fan.RAYS))), dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, scalars)
    stacks = Dict(pairs)

    return(StackyFan(fan, stacks))
end

"""
    addStackStructure(::Polymake.BigObjectAllocated, ::Array{Int64, 1})

Function to generate a stacky fan from a given fan and a set of scalars.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 1; 1 2],INPUT_CONES=[[0,1],[1,2]]);

julia> stackyWeights(addStackStructure(X,[2,2,2]))
[ 2 ,  2 ,  2 ]
```
"""
function addStackStructure(
    fan::Polymake.BigObjectAllocated,
    scalars::Array{Int64, 1})
    
    # Construct the dictionary
    stack_rays = mapslices(encode, Polymake.common.primitive(fan.RAYS), dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, scalars)
    stacks = Dict(pairs)

    return(StackyFan(fan, stacks))
end

"""
    encode(::Union{Polymake.VectorAllocated{Polymake.Rational},Polymake.VectorAllocated{Polymake.Integer},Vector{Int64}})

    Internal function that converts a vector, representing a ray in the fan,
to a string in order to allow for hashing for the dictionary.

#Examples
```jldoctest
julia> encode([1,0,2,5])
"1,0,2,5"
```
"""
function encode(objects::Union{Polymake.VectorAllocated{Polymake.Rational},Polymake.VectorAllocated{Polymake.Integer},Vector{Int64}})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

"""
    stackyWeights(::StackyFan)

    Returns a list of the stacky weights of the rays of the given stacky fan with the same order as the rays of the fan.

#Examples
```jldoctest
julia> F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,2,3,4]);

julia> stackyWeights(F)
[ 1 ,  2 ,  3 ,  4 ]
```
"""
function stackyWeights(sf::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(sf.fan.RAYS)))
    #println("Stacky weights ray matrix: $rayMatrix")
    rayList=slicematrix(rayMatrix)
    out=Int64[]
    for ray in rayList
        stack=sf.stacks[encode(ray)]
        push!(out,stack)
    end
    return out
end

## API functions

"""
    getRayStack(::StackyFan, ::Array{Int64, 1})

    Get the scalar associated with a ray in the given stacky fan structure.

# Examples
```jldoctest
julia> F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,2,3,4]);

julia> getRayStack(F,[1,2])
3
```
"""
function getRayStack(sf::StackyFan, ray::Array{Int64, 1})
    return sf.stacks[encode(ray)]
end

"""
    rootConstruction(::StackyFan, ::Array{Int64, 1})

Given a fan and a set of scalars corresponding to the rays of the fan,
performs a root construction on the fan by multiplying the stack scalars
by the given values. 

rootConstruction returns a new StackyFan object, and does not modify the input.

# Examples
```jldoctest

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);

julia> SX = StackyFan(X, [2,3,5,7]);

julia> stackyWeights(rootConstruction(SX, [1, 4, 2, 1]))
[ 2 ,  12 ,  10 ,  7 ]
```
"""
function rootConstruction(
    sf::StackyFan,
    scalars::Array{Int64, 1})
    
    # Multiply the scalars of the fan by the given values
    return StackyFan(sf.fan, stackyWeights(sf) .* scalars)
end

"""
    rootConstructionDistinguishedIndices(::StackyFan, ::Array{Int64, 1}, ::Array{Int64, 1})

    Given a fan, the indices of the distinguished rays in the fan rays (as an incidence matrix), and
a set of scalars corresponding to the rays of the fan, performs a root 
construction on the fan by multiplying the stack scalars by the given values. 

    rootConstructionDistinguishedIndices returns a new StackyFan object,
and does not modify the input.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);

julia> SX = StackyFan(X, [2,3,5,7]);

julia> stackyWeights(rootConstructionDistinguishedIndices(SX, [0,1,1,0], [4, 2, 1, 3]))
[ 2 ,  6 ,  5 ,  7 ]
```
"""
function rootConstructionDistinguishedIndices(
    sf::StackyFan,
    distIndices::Array{Int64, 1},
    scalars::Array{Int64, 1})
    
    numRays = size(sf.fan.RAYS, 1)
    fullScalars = fill(1, numRays)
    for i in 1:numRays
        if distIndices[i]==1 && scalars[i] != 0
            fullScalars[i] = scalars[i]
        end
    end
    # Multiply the scalars of the fan by the given values
    return rootConstruction(sf, fullScalars)
end

"""
    rootConstructionDistinguished(
        ::StackyFan, 
        ::Polymake.Matrix{Polymake.Rational},
        ::Array{Int64, 1})

    Given a fan, a set of distinguished rays, and a set of scalars of equal size,
performs a root construction on the fan on the distinguished rays by multiplying 
the stack scalars by the given values.

    rootConstructionDistinguished returns a new StackyFan object, 
and does not modify the input.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);

julia> SX = StackyFan(X, [2,3,5,7]);

julia> distinguished = X.RAYS[[2,3],:];

julia> stackyWeights(rootConstructionDistinguished(SX, distinguished, [4, 2]))
[ 2 ,  12 ,  10 ,  7 ]
```
"""
function rootConstructionDistinguished(
    sf::StackyFan,
    rays::Polymake.Matrix{Polymake.Rational},
    scalars::Array{Int64, 1})
    
    # Check that the rays and scalars are the same size
    #if (size(rays, 1) != length(scalars))
    #    error("Inputs are not of equal size")
    #end
    
    encoded_rays = mapslices(encode, rays, dims=2)
    # Make a copy of the dictionary
    newStacks = copy(sf.stacks)
    for i in 1:length(encoded_rays)
        ray = encoded_rays[i]
        # Multiply the scalar of the corresponding ray
        newStacks[ray] *= scalars[i]
    end
    
    # Convert the dictionary to an array of scalars matching the indices
    #newScalars = mapslices(ray -> newStacks[encode(ray)], sf.fan.RAYS, dims=2)
    newScalars = Array{Int64, 1}()
    for i in 1:size(sf.fan.RAYS, 1)
        push!(newScalars, newStacks[encode(sf.fan.RAYS[i,:])])
    end
    
    return StackyFan(sf.fan, newScalars)
end

"""
    findBarycenter(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated)

    Takes a normal toric variety X and a set s corresponding to a subset of rays of X, and outputs a polymake vector corresponding to the barycenter of those rays.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);

julia> s=[1,2];

julia> findBarycenter(s,X)
pm::Matrix<pm::Integer>
2 1
```
"""
function findBarycenter(s::Union{AbstractSet,AbstractVector},X::Polymake.BigObjectAllocated)
#     rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
#     rays = rowMinors(rayMatrix, s)
#     dim=size(rays,2)
#     bary=zeros(Int64,dim,1)
#     for i in 1:size(rays,1)
#         bary+=rays[i,:]
#     end
#     bary=Polymake.common.primitive(bary)
#     return vec(bary)
# end
    
    rays = rowMinors(Array(X.RAYS), s)
    dim=size(rays,2)
    bary=zeros(Polymake.Rational,dim,1)
    for i in 1:size(rays,1)
        bary+=rays[i,:]
    end
    bary=Polymake.common.primitive(transpose(bary))
    return bary
end

"""
    findStackyBarycenter(::Union{AbstractSet,AbstractVector},::StackyFan)

    Takes a stacky fan SX and a set s corresponding to a subset of rays of SX, calculates the 'stacky rays' corresponding to those rays (the rays times their stacky values), and find the barycenter of the stacky rays.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[2,3]);

julia> findStackyBarycenter([1,2],F)
[ 5 ,  6 ]
```
"""
function findStackyBarycenter(s::Union{AbstractSet,AbstractVector},SX::StackyFan)
    rayMatrix=convert(Array{Int64,2}, Array(Polymake.common.primitive(SX.fan.RAYS)))
    # Multiply the rays by their stacky values
    stackMatrix = diagm(stackyWeights(SX)) * rayMatrix
    rays = rowMinors(stackMatrix, s)
    dim=size(rays,2)
    bary=zeros(Int64,dim,1)
    for i in 1:size(rays,1)
        bary+=rays[i,:]
    end
    return vec(bary)
end

"""
    findStackyRayMatrix(::StackyFan)

    Outputs the ray matrix of the given stacky fan, such that each primitive ray is multiplied by its stacky weights.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 3; 1 5],INPUT_CONES=[[0,1],[1,2]]);

julia> F=addStackStructure(X,[1,2,3]);

julia> findStackyRayMatrix(F)
3×2 Matrix{Int64}:
 1   0
 4   6
 3  15
```
"""
function findStackyRayMatrix(sf::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(sf.fan.RAYS)))
    sW=stackyWeights(sf)
    return sW .* rayMatrix
end


"""
    getConeRank(::AbstractMatrix,::AbstractVector)

Takes a matrix and a vector containing indices corresponding to rows of a matrix, and calculates the rank of the matrix consisting only of those rows.

#Examples
```jldoctest makeSmoothWithDependencies
julia> v=[1,2]

julia> M=[0 1; 1 1; 1 0]

julia> getConeRank(v,M)
2
```
"""
function getConeRank(coneRayIndices::AbstractVector, rayMatrix::AbstractMatrix)
    coneRays = rowMinors(rayMatrix,coneRayIndices)
    return rank(Matrix(coneRays))
end

"""
    getDimension(::Polymake.BigObjectAllocated)

Returns the ambient dimension of a normal toric variety.

#Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 2 0;0 0 1;0 1 0; 1 1 1],INPUT_CONES=[[0,1,2],[0,2,3,4]])

julia> getDimension(X)
3
```
"""
function getDimension(X)
    return size(X.RAYS, 2)
end

"""
    getConeFaces(::Polymake.BigObjectAllocated,::AbstractVector,::AbstractMatrix)

Takes a fan, its ray matrix, and a vector corresponding to one of its cones, and returns a list of maximal strict faces of that cone.

#Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 1],INPUT_CONES=[[0,1,2,3]])

julia> getConeFaces(X,[1,2,3,4],Array(X.RAYS))
[[ 1 ,  2 ], [ 1 ,  3 ], [ 3 ,  4 ], [ 2 ,  4 ]]
```
"""
function getConeFaces(fan::Polymake.BigObjectAllocated,cone::AbstractVector,rayMatrix::AbstractMatrix)
    lattice = fan.HASSE_DIAGRAM
    faces = @Polymake.convert_to Array{Set{Int}} lattice.FACES
    cone_faces=[]
    c = rank(Array(rowMinors(rayMatrix, cone))) - 1
    rank_c_subcone_indices = @Polymake.convert_to Array{Int} Polymake.graph.nodes_of_rank(lattice,c)
    rank_c_subcones = [faces[i + 1] for i in rank_c_subcone_indices]
    for subcone in rank_c_subcones
        new_cone = [i+1 for i in subcone]
        if all((i -> i in cone).(new_cone))
            push!(cone_faces, new_cone)
        end
    end 
    return cone_faces
end


"""
    toric_blowup(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated,::AbstractVector)

    Takes a normal toric variety X, a set s corresponding to a subset of rays of X, and a (optional) polymake vector, v, blow up X at v. If v is not provided, blow up X at the barycenter of s.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);

julia> B=toric_blowup([0,1],X,[1,1]);

julia> B.INPUT_RAYS
pm::Matrix<pm::Rational>
1 0
1 2
1 1

julia> B.INPUT_CONES
pm::IncidenceMatrix<pm::NonSymmetric>
{0 2}
{1 2}
```
"""
function toric_blowup(s, X, v)
    s = [i + 1 for i in s]
    if v==nothing
        v=findBarycenter(s,X)
    end
    if size(v,2)==1
         v=transpose(v)
    end
    coneList = convertIncidenceMatrix(X.MAXIMAL_CONES)
    # Extracting the indices of all the cones in X that contains the set of rays s
    starIndex = findall((t) -> all(((i) -> i in t).(s)), coneList)
    star = [coneList[i] for i in starIndex]
    rayMatrix = X.RAYS
    
    lattice = X.HASSE_DIAGRAM
    faces = @Polymake.convert_to Array{Set{Int}} lattice.FACES
    
    # Get all the subcones of X that is contained in one of the cones in star and has the same rank
    clStar = []
    # Iterate over star
    for t in star
        # Get the rank of t
        c = rank(Array(rowMinors(rayMatrix, t))) - 1
        # Get all the subcones of X with rank c
        rank_c_subcone_indices = @Polymake.convert_to Array{Int} Polymake.graph.nodes_of_rank(lattice,c)
        rank_c_subcones = [faces[i + 1] for i in rank_c_subcone_indices]
        # Iterate over rank_c_subcones, and put the cones that is contained in t into clStar
        for cone in rank_c_subcones
            new_cone = [i+1 for i in cone]
            if all((i -> i in t).(new_cone))
                push!(clStar, new_cone)
            end
        end
    end
    clStar = unique(clStar)
    
    n = size(rayMatrix, 1) + 1
    # Filter out the cones in star from conelist
    coneList = filter(x -> !(x in star), coneList)
    
    if length(s) == 1
        # If s consists of a single ray, find all the cones in clStar that does not contain s
        newCones = []
        for t in clStar
            if !(s[1] in t)
                push!(newCones, sort(push!(t, s[1])))
            end
        end
        # return newCones plus coneList
        finalCones = [[i - 1 for i in cone] for cone in append!(coneList, newCones)]
        return Polymake.fulton.NormalToricVariety(INPUT_RAYS = Array(X.RAYS), INPUT_CONES = finalCones)
    end
    newCones = []
    for t in clStar
        # Find all the cones in clStar that does not contain at least one ray in s
        # QUESTION: Why seperate this from the one element case? Any won't work with one element list?
        if any(((i) -> !(i in t)).(s))
            push!(newCones, push!(t, n))
        end
    end
    # return newCones plus coneList
    finalRays = vcat((X.RAYS),v)
    finalCones = [[i - 1 for i in cone] for cone in append!(coneList, newCones)]
    return Polymake.fulton.NormalToricVariety(INPUT_RAYS = finalRays, INPUT_CONES = finalCones)
end

"""
    makeSimplicial(::Polymake.BigObjectAllocated)

Takes in a normal toric variety and returns a simplicial toric variety by subdividing (blowing up) the non-simplicial maximal cones.

#Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])

julia> X.SIMPLICIAL
false

julia> makeSimplicial(X).SIMPLICIAL
true

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0;0 1 0 0;0 0 1 0;1 -1 1 0; 1 0 -2 0],INPUT_CONES=[[0,1,2,3],[0,4]])

julia> X.SIMPLICIAL
false

julia> makeSimplicial(X).SIMPLICIAL
true
```
"""
function makeSimplicial(X::Polymake.BigObjectAllocated)
    Y = copy(X)
    while (true)
        # If the initial toric variety is simplicial, the program terminates and returns it.
        if Y.SIMPLICIAL==true
            break
        end
        #Maximal cones and ray matrix
        coneList = convertIncidenceMatrix(Y.MAXIMAL_CONES)
        rayMatrix = Y.RAYS
        badCone = nothing
        for i in 1:size(coneList,1)
            cone = coneList[i]
            if (getConeRank(cone, rayMatrix) != size(cone)[1])
                badCone = cone
            end
        end
        if (badCone == nothing)
            # All cones are linearly independent
            break
        else
            # Find the first ray that is contained in more than one orbit
            # and subdivide at that ray, using toricBlowup
            
            # Get faces (need to replace this)
            edges = getConeFaces(Y,badCone,rayMatrix)
            # Find the first ray that is contained in more than one orbit
            i = 1
            while count(r->(badCone[i] in r), edges) == 1
                i += 1
            end
            # Subdivide at the cone containing just that ray
            badCone=[i-1 for i in badCone]
            Y = toric_blowup(badCone, Y,nothing)
            #Y = toric_blowup([badCone[i]], Y,nothing)
        end
        # Repeat this process until there are no more bad cones
    end
    return Y
end

"""
    makeSmooth(::Polymake.BigObjectAllocated)

Takes in a normal toric variety X and output a new smooth toric variety by iteratively blowing up. In the language of fans, these blowups are achieved by subdividing non-smooth cones.

#Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[4 -1; 0 1],INPUT_CONES=[[0, 1]])

julia> X.SMOOTH_FAN
false

julia> makeSmooth(X).SMOOTH_FAN
true

julia> makeSmooth(X).INPUT_RAYS
pm::Matrix<pm::Rational>
1 -1/4
0 1
1 0

julia>X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])

julia>X.SMOOTH_FAN
false

julia>makeSmooth(X).SMOOTH_FAN
true

julia>makeSmooth(X).INPUT_RAYS
pm::Matrix<pm::Rational>
1 0 0
1 1 0
1 0 1
1 1 1
2 1 1

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;0 1 0;0 0 1;0 -1 -1; -1 0 -1; -2 -1 0],INPUT_CONES=[[0,1,2],[0,1,3],[1,3,4],[1,2,4],[2,4,5],[0,2,5],[0,3,5],[3,4,5]])

julia> X.SMOOTH_FAN
false

julia> makeSmooth(X).SMOOTH_FAN
true
```
"""
function makeSmooth(X::Polymake.BigObjectAllocated)
    Y  = copy(X)
    while(true)
        coneList = convertIncidenceMatrix(Y.MAXIMAL_CONES)
        rayMatrix = Array(Polymake.common.primitive(Y.RAYS))
        k = 1
        # Iterate through the coneList, getting the index of the first cone not smooth
        for coneSet in coneList
            # Getting the number of rays in coneSet
            S=size(coneSet)[1]
            coneRays=rowMinors(rayMatrix,coneSet)
            # Checking whether this cone is smooth
            smoothCheck=Polymake.fan.check_fan_objects(Polymake.polytope.Cone(RAYS=coneRays)).SMOOTH_FAN
            if !smoothCheck
                # If the cone is not simplicial or not smooth, we have found the cone that we need to make smooth
                break
            else
                k+=1
            end
        end
        # At this point, all the cones are smooth. The program terminates.
        if k == size(coneList,1)+1
            break
        end
        
        # Get the cone that we found to be not smooth
        sigma=coneList[k]
        sigmaRays=slicematrix(rowMinors(rayMatrix,sigma))
        tau=0; tauRays=0; tauCone=0
        # Iterate over the subcones of sigma, finding tau, the smallest one that is not smooth
        for subset in collect(powerset(sigma))
            if size(subset,1) > 1
                S=size(subset)[1]
                subsetRays=rowMinors(rayMatrix,subset)
                subsetCone=Polymake.polytope.Cone(RAYS=subsetRays)
                smoothCheck=Polymake.fan.check_fan_objects(subsetCone).SMOOTH_FAN
                if !smoothCheck
                    tau=subset
                    tauRays=subsetRays
                    tauCone=subsetCone
                    break
                end 
            end
        end
        
        # Getting the Hilbert Basis of tau
        H=slicematrix(Matrix(tauCone.HILBERT_BASIS_GENERATORS[1]))
        rayIndex=0
        # Iterate over the Hilbert Basis, finding the first ray that is not the generator of sigma
        for i in 1:size(H,1)
            if !(H[i] in sigmaRays)
                rayIndex=i
                break
            end
        end
        if rayIndex==0
            # Every Hilbert Basis of tau is a generator of sigma. Make Y simplicial is sufficient to make sigma smooth
            Y=makeSimplicial(Y)
        else
            # blowupRay is not a generator of sigma, blow up tau at blowupRay
            blowupRay=H[rayIndex]
            tau=[i-1 for i in tau]
            Y=toric_blowup(tau,Y,transpose(blowupRay))
        end
    end
    return Y
end

"""
    stackyBlowup(::StackyFan,::Array{Int64,1},::Array{Int64,1})

    Takes a stacky fan sf, a ray excep, and a cone, and subdivides the stacky fan at the given ray. Crucially, the given cone should be the minimal cone containing the exceptional ray. The cone input should be zero-indexed.

#examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[2,3]);

julia> stackyWeights(stackyBlowup(F,[0,1],[1,1]))
[ 2 ,  1 ,  3 ]
```
"""
function stackyBlowup(sf::StackyFan, cone::Array{Int64,1}, excep::Array{Int64,1})
    # Express the exceptional ray as a scalar multiple of a primitive ray
    # Use this scalar as the stacky weight in the resulting stacky fan
    G=gcd(excep)
    excep=Polymake.common.primitive(excep)
    
    # Perform toric blowup at the given ray
    blowup = toric_blowup(cone, sf.fan, excep)
    sf.stacks[encode(excep)] = G

    return(StackyFan(blowup, sf.stacks))
end

"""
    getConesPolymake(sf::StackyFan)

    Returns a list of cones of a stacky fan, with the cones represented as polymake objects.
    
#Examples
```jldoctest
julia> F=makeStackyFan([1 0; 1 1; 1 2; 1 3],[[0,1],[1,2],[2,3]],[1,1,1,1])
    
julia> getConesPolymake(F)[1].RAYS
pm::Matrix<pm::Rational>
1 0
1 1
```
"""
function getConesPolymake(sf::StackyFan)
    formatted = convertIncidenceMatrix(sf.fan.MAXIMAL_CONES)
    cones = map((x) -> Polymake.polytope.Cone(
        INPUT_RAYS=sf.fan.RAYS[x,:]), formatted)
    return(cones)
end


"""
    slicematrix(::AbstractMatrix{<:Number})

    Take a two-dimensional matrix and output a list of its row vectors.

# Examples
```jldoctest
julia> A=[1 2; 3 4];

julia> slicematrix(A)
[[ 1 ,  2 ], [ 3 ,  4 ]]
```
"""
function slicematrix(A::AbstractMatrix{<:Number})
    return [A[i, :] for i in 1:size(A,1)]
end

"""
    rowMinors(::AbstractMatrix{<:Number},::Union{AbstractSet,AbstractVector})

    Identical to slicematrix, except only returns row vectors indexed by a set S.

# Examples
```jldoctest
julia> A=[1 2 3;4 5 6; 7 8 9];

julia> S=Set([1,3]);

julia> rowMinors(A,S)
2×3 LinearAlgebra.Transpose{Int64,Array{Int64,2}}:
 1  2  3
 7  8  9
```
"""
function rowMinors(A::AbstractMatrix{<:Number},S::Union{AbstractSet,AbstractVector})
    outList=[]
    slices=slicematrix(A)
    for i in 1:size(slices,1)
        if i in S
            append!(outList,[slices[i]])
        end
    end
    return Array(transpose(hcat(outList...)))
end

"""
    convertIncidenceMatrix(::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})

    Takes a Polymake incidence matrix (e.g., the output of X.MAXIMAL_CONES for a toric variety X) and outputs a list of vectors, with each vector recording the indices marked on a given row of the incidence matrix.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);

julia> M=X.MAXIMAL_CONES;

julia> convertIncidenceMatrix(M)
[[ 1 ,  2 ], [ 2 ,  3 ]]
```
"""
function convertIncidenceMatrix(A::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})
    A=Array(A)
    dim1=size(A,1)
    dim2=size(A,2)
    out=[]
    for i in 1:dim1
        members=[]
        for j in 1:dim2
            if A[i,j]==true
                append!(members,j)
            end
        end
        append!(out,[members])
    end
    return convert(Array{Array{Int64,1},1}, out)
end

"""
    coneMultiplicity(C::Polymake.BigObjectAllocated)

    Returns the multiplicity of a polyhedral cone (inputted as a Polymake object): here, the multiplicity is defined as the index of the sublattice generated by the rays of the cone, inside the full integer lattice contained in the linear subspace generated by the edges of the cone.

# Examples
```jldoctest

julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 0; 1 2]);

julia> coneMultiplicity(C)
2
```
"""
function coneMultiplicity(C::Polymake.BigObjectAllocated)
    A=Polymake.common.primitive(C.RAYS)
    M=matrix(ZZ,[fmpz.(y) for y in A])
    SNF=Nemo.snf(M)
    mult=1
    for i in 1:size(SNF,1)
        mult*=SNF[i,i]
    end
    return mult
end

"""
    getMultiplicity(::Array{Int64,1},::Array{Int64,2})
        
    Same functionality as coneMultiplicity, but calculates the cone rays as a subset of the columns of a ray matrix rather than from a Polymake cone object.
        
# Examples
```jldoctest
julia> getMultiplicity([1,2],[1 0; 1 2; 1 3])
2       
```     
"""
function getMultiplicity(cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    A=rowMinors(rayMatrix,cone)
    M=matrix(ZZ,[fmpz.(y) for y in A])
    SNF=Nemo.snf(M)
    mult=1
    for i in 1:size(SNF,1)
        mult*=SNF[i,i]
    end
    return mult
end

"""
    coneConvert(::abstractVector{Int64},::abstractMatrix{Int64})

    Takes a matrix where the columns represent rays, and a list of indices, and forms a Polymake cone object generated by the rays corresponding to those indices.

# Examples
```jldoctest

julia> typeof(coneConvert([1, 2, 4],[1 0 0; 0 1 0; 0 0 1; 1 1 1]))
Polymake.BigObjectAllocated
```
"""
function coneConvert(cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    coneRays=rowMinors(rayMatrix,cone)
    C=Polymake.polytope.Cone(RAYS=coneRays)
    return C
end
    
"""
    getCones(X::Polymake.BigObjectAllocated)
    
    Returns all the cones of a fan X as a list of lists, with each interior list containing the indices of the rays generating a given cone.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 ],INPUT_CONES=[[0,1,2,3]]);

julia> getCones(X)
[[ 0 ,  1 ,  2 ,  3 ], [ 0 ,  1 ], [ 0 ,  2 ], [ 2 ,  3 ], [ 1 ,  3 ], [ 0 ], [ 1 ], [ 2 ], [ 3 ]]
```
"""
function getCones(X::Polymake.BigObjectAllocated)
    lattice=X.HASSE_DIAGRAM
    faces=@Polymake.convert_to Array{Set{Int}} lattice.FACES
    out=Array{Int64,1}[]
    for i in 2:(size(faces,1)-1)
        newface=Array(@Polymake.convert_to Array{Int} faces[i])
        push!(out,[i+1 for i in newface])
    end
    return out
end

"""
    findFaceContainingRay(::Polymake.BigObjectAllocated,::Array{Int64, 1})

Given a cone and a ray, finds a maximal proper face of the cone containing that ray.

# Examples
```jldoctest
julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1]);

julia> findFaceContainingRay(C,[1,1,1]).RAYS
pm::Matrix<pm::Rational>
1 0 1
1 1 1
```
"""
function findFaceContainingRay(C::Polymake.BigObjectAllocated,v::Array{Int64, 1})
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(C.RAYS)))
    faces=convertIncidenceMatrix(C.RAYS_IN_FACETS)
    if faces==[[]]
        return nothing
    end
    for face in faces
        faceCone=coneConvert(face,rayMatrix)
        if Polymake.polytope.contains(faceCone, v)
            return faceCone
        end
    end
    return nothing
end

"""
    findMinimalCone(::Polymake.BigObjectAllocated,::Array{Int64,1})

Given a fan and a ray, finds the minimal cone fo the fan containing the ray.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0; 1 0 1 0; 1 1 0 0; 1 1 1 0; 0 0 0 1],INPUT_CONES=[[0,1,2,3,4]]);

julia> findMinimalCone(X,[1,2,0,0])
[1, 3]
```
"""
function findMinimalCone(X::Polymake.BigObjectAllocated,v::Array{Int64, 1})
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    cones=convertIncidenceMatrix(X.MAXIMAL_CONES)
    #print(cones)
    startCone=nothing
    for cone in cones
        #print(cone)
        polyCone=coneConvert(cone,rayMatrix)
        if Polymake.polytope.contains(polyCone, v)
            startCone=polyCone
        end
    end
    if startCone==nothing
        error("The given ray is not contained in any cone of the fan.")
    end
    currentCone=startCone
    while(true)
        nextCone=findFaceContainingRay(currentCone,v)
        if nextCone==nothing
            break
        end
        currentCone=nextCone
    end
    currentRays=slicematrix(Array(currentCone.RAYS))
    fanRays=slicematrix(rayMatrix)
    indices=findall(x->x in currentRays,fanRays)
    return indices
end

"""
    starSubdivision(::Polymake.BigObjectAllocated,::Array{Int64,1})

    Blows up the input fan at the input ray. A wrapper integrating findMinimalCone and toric_blowup.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0 0; 1 0 1 0; 1 1 0 0; 1 1 1 0; 0 0 0 1],INPUT_CONES=[[0,1,2,3,4]]);

julia> starSubdivision(X,[2,1,0,0]).RAYS
pm::Matrix<pm::Rational>
1 0 0 0
1 0 1 0
0 0 0 1
1 1/2 0 0
1 1 0 0
1 1 1 0
```
"""
function starSubdivision(X::Polymake.BigObjectAllocated, v::Array{Int64, 1})
    minimalCone = findMinimalCone(X, v)
    s = [i - 1 for i in minimalCone]
    v = transpose(v)
    return toric_blowup(s, X, v)
end

"""
    convertToIncidence(v::Array{Int64,1},l::Int64)

Returns a vector of length l, with entries of 1 indexed by v and entries of 0 everywhere else.

# Examples
```jldoctest
julia> convertToIncidence([2,3,5],6)
[ 0 , 1 , 1 , 0 , 1 , 0 ]
```
"""
function convertToIncidence(v::Array{Int64,1},l::Int64)
    out=[]
    for j in 1:l
        if j in v
            append!(out,1)
        else
            append!(out,0)
        end
    end
    return out
end

"""
    interiorPoints(::Polymake.BigObjectAllocated)

    Finds all interior lattice points contained in the fundamental region of a given cone. When multiple interior lattice points lie along the same ray, only the point closest to the origin is returned. Notably, 

# Examples
```jldoctest
julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 2; 2 1]);

julia> interiorPoints(C)
[[ 1 ,  1 ]]
```
"""
function interiorPoints(C::Polymake.BigObjectAllocated)
    rayMatrix=Array(Polymake.common.primitive(C.RAYS))
    l=size(rayMatrix,1)
    dim=size(rayMatrix,2)
    if rank(rayMatrix)<l
        print(rayMatrix)
        error("Input cone is not simplicial.")
    end
    subsets=collect(powerset([1:l;]))
    vertices=[]
    for elt in subsets #vertices of the fundamental region are in correspondence with subsets of the generators of the cone, by summing the generators in a subset to obtain a vertex
        vert=zeros(Polymake.Rational,1,dim)
        for i in 1:l
            if i in elt
                vert+=rayMatrix[[i],:]
            end
        end
        append!(vertices,[vert])
    end
    V=vcat(vertices...)
    VH=hcat(ones(Polymake.Rational,size(V,1)),V)
    P=Polymake.polytope.Polytope(POINTS=VH) #make a Polymake polytope object from the vertices of the fundamental region found in the last step
    if size(P.INTERIOR_LATTICE_POINTS,1)==0
        return nothing
    end
    intPoints=Array(P.INTERIOR_LATTICE_POINTS)[:,2:(dim+1)] #find all the interior lattice points
    validPoints=[]
    #return intPoints
    for i in 1:size(intPoints,1) #throw out all points that are integer multiples of other points
        point=intPoints[i,:]
        if gcd(point)==1
            append!(validPoints,[point])
        end
    end
    return validPoints
end

"""
    coneRayDecomposition(::Array{Int64,1},::Array{Int64,2},::Array{Int64,1},::Array{Int64,1})

    This function takes in a cone (a vector of indices of cone generators in rayMatrix), a ray, and a stacky structure for rayMatrix. It first multiplies all generators of the cone by their stacky values, and then finds an expression for the ray as a sum of these stacky generators. The output is a vector of coefficients of the above representation in terms of the rays in rayMatrix, with zeros as coefficients for all rays not in the given cone.

# Examples
```jldoctest
julia> coneRayDecomposition([1,2,3],[3 5 7; 8 16 9;2 1 3;1 1 1],[2,2,3],[1,1,1,1])
[ 6 ,  5 ,  52 ,  0 ]
```
"""
function coneRayDecomposition(cone,rayMatrix,ray,stack)
    stackMatrix=diagm(stack)*rayMatrix # multiply all rays by stack values
    coneRays=rowMinors(stackMatrix,cone) # find the (stackified) rays in the given cone
    if rank(coneRays)<size(coneRays,1)
        error("The given cone is not simplicial.")
    end
    B=Polymake.common.null_space(hcat(transpose(coneRays),-ray)) # Express the input ray in terms of the stackified cone generators
    N=convert(Array{Int64,1},vec(B))
    if size(N,1)==0
        error("The given ray is not in the span of the cone generators.")
    end
    if N[end]<0 #since the nullspace has arbitrary sign, fix it so the coefficients are all positive
        N*=-1
    end
    pop!(N)
    out=zeros(Int64,size(rayMatrix,1)) 
    for i in 1:size(N,1)#rewrite the coefficients vector in terms of all the rays in rayMatrix, by padding with zeros when appropriate.
        out[cone[i]]=N[i] 
    end
    return out
end

"""
=========VISUALIZATION FUNCTIONALITY==========
"""

using Plots

"""
    plot3dSimpCone(::Array{Array{Int64,1},1})

    Give a list of three 3-dimensional vectors, plots the polygon defined by those vectors and the origin via Plots.jl.

# Examples
```jldoctest
julia> plot3dSimpCone([[1,0,0],[0,1,0],[0,0,1]]);
```
"""
function plot3dSimpCone(L::Array{Array{Int64,1},1})
    l=size(L,1)
    if l==3
        (x1,y1,z1)=L[1]
        (x2,y2,z2)=L[2]
        (x3,y3,z3)=L[3]
        X=[0,x1,0,x2,0,x3,x2,x1,x3]
        Y=[0,y1,0,y2,0,y3,y2,y1,y3]
        Z=[0,z1,0,z2,0,z3,z2,z1,z3]
        plot!(X,Y,Z)
    elseif l==2
        (x1,y1,z1)=L[1]
        (x2,y2,z2)=L[2]
        X=[0,x1,x2,0]
        Y=[0,y1,y2,0]
        Z=[0,z1,z2,0]
        plot!(X,Y,Z)
    elseif l==1
        (x1,y1,z1)=L[1]
        X=[0,x1]
        Y=[0,y1]
        Z=[0,z1]
        plot!(X,Y,Z)
    else 
        error("Input cone is empty or not simplicial.")
    end
end

"""
    plot2dCone(::Array{Array{Int64,1},1})

    Given a list of two 2-dimesnional vectors, plots the polygon defined by those vectors and the origin via Plots.jl.

# Examples
```jldoctest
julia> plot2dCone([[1,0],[0,1]]);
```
"""
function plot2dCone(L::Array{Array{Int64,1},1})
    l=size(L,1)
    if l==2
        (x1,y1)=L[1]
        (x2,y2)=L[2]
        S=Shape([(0,0),(x1,y1),(x2,y2)])
        plot!(S)
    elseif l==1
        (x1,y1)=L[1]
        X=[0,x1]
        Y=[0,y1]
        plot!(X,Y)
    else
        error("Input cone is empty or not simplicial.")
    end
end

"""
    showSimpFan(X::Polymake.BigObjectAllocated)

    Plots the input simplicial fan as a collection of polygons or polyhedra defined by the input rays and the origin via Plots.jl.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);

julia> showSimpFan(X);
```
"""
function showSimpFan(X::Polymake.BigObjectAllocated)
    plot()
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    dim=size(rayMatrix,2)
    if !(dim in [2,3])
         println("Fan visualization is only supported in dimensions 2 and 3.")
    end
    maxCones=convertIncidenceMatrix(X.MAXIMAL_CONES)
    for cone in maxCones
        coneRays=slicematrix(rowMinors(rayMatrix,cone))
        if dim==2
            plot2dCone(coneRays)
        elseif dim==3
            plot3dSimpCone(coneRays)
        end
    end
    plot!(legend=false)
end

"""
    showSimpStackyFan(::StackyFan;::Bool=true)

    Extends the functionality of showSimpFan to stacky simplicial fans. If the stackypoints input is set to true, a red dot is shown at the location of each primitive vector multiplied by its stacky weight.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,2]);

julia> showSimpStackyFan(F);
```
"""
function showSimpStackyFan(F::StackyFan;stackypoints::Bool=true)
    X=F.fan
    plot()
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    dim=size(rayMatrix,2)
    if !(dim in [2,3])
         println("Fan visualization is only supported in dimensions 2 and 3.")
    end
    maxCones=convertIncidenceMatrix(X.MAXIMAL_CONES)
    for cone in maxCones
        coneRays=slicematrix(rowMinors(rayMatrix,cone))
        if dim==2
            plot2dCone(coneRays)
        elseif dim==3
            plot3dSimpCone(coneRays)
        end
    end
    if stackypoints==true
        stackyPoints=slicematrix(findStackyRayMatrix(F))
        for point in stackyPoints
            if dim==2
                (x,y)=point
                plot!([x],[y],shape=:circle,color=:red)
            elseif dim==3
                (x,y,z)=point
                plot!([x],[y],[z],shape=:circle,color=:red)
            end
        end
    end
    plot!(legend=false)
end

"""
    coneVectorOrder(::Polymake.BigObjectAllocated)

    Takes a polyhedral cone in 2 or 3 dimensions and outputs a list of its defining rays, arranged in counterclockwise order around the exterior of the cone.

# Examples
```jldoctest
julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1]);

julia> coneVectorOrder(C)
[[1
, 1, 0], [1, 1, 1], [1, 0, 1], [1, 0, 0]]
```
"""
function coneVectorOrder(C::Polymake.BigObjectAllocated)
    rayList=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(C.RAYS))))
    orderedRayList=Array{Int64,1}[]
    ordering=Array(C.RIF_CYCLIC_NORMAL[1])
    ordering=[i+1 for i in ordering]
    for i in ordering
        push!(orderedRayList,rayList[i])
    end
    return orderedRayList
end

"""
    plot3dCone(::Array{Array{Int64,1},1})

    Plots the 3-dimensional cone defined by the given list of rays. The input rays are assumed to be in counterclockwise or clockwise order around the cone; coneVectorOrder() can be used to obtain this ordering from an arbitrary cone.

# Examples
```jldoctest
julia> plot3dCone([[1,1,0],[1,1,1],[1,0,1],[1,0,0]]);
```
"""
function plot3dCone(L::Array{Array{Int64,1},1})
    l=size(L,1)
    X=Int64[0]
    Y=Int64[0]
    Z=Int64[0]
    for i in 1:l
        push!(X,L[i][1])
        push!(Y,L[i][2])
        push!(Z,L[i][3])
        push!(X,0)
        push!(Y,0)
        push!(Z,0)
    end
    for i in 1:l
        push!(X,L[i][1])
        push!(Y,L[i][2])
        push!(Z,L[i][3])
    end
    push!(X,L[1][1])
    push!(Y,L[1][2])
    push!(Z,L[1][3])
    plot!(X,Y,Z)
end

"""
    showFan(X::Polymake.BigObjectAllocated)

    Plots the input fan as a collection of polygons or polyhedra defined by the input rays and the origin via Plots.jl.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1],INPUT_CONES=[[0,1,2,3]]);

julia> showFan(X);
```
"""
function showFan(X::Polymake.BigObjectAllocated)
    plot()
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    dim=size(rayMatrix,2)
    if !(dim in [2,3])
         println("Fan visualization is only supported in dimensions 2 and 3.")
    end
    maxCones=convertIncidenceMatrix(X.MAXIMAL_CONES)
    for cone in maxCones
        if dim==2
            coneRays=slicematrix(rowMinors(rayMatrix,cone))
            plot2dCone(coneRays)
        elseif dim==3
            C=coneConvert(cone,rayMatrix)
            plot3dCone(coneVectorOrder(C))
        end
    end
    plot!(legend=false)
end
    
"""
    showStackyFan(::StackyFan;::Bool=true)

    Extends the functionality of showFan to stacky fans. If the stackypoints input is set to true, a red dot is shown at the location of each primitive vector multiplied by its stacky weight.

# Examples
```jldoctest
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 0 1; 1 1 0; 1 1 1],INPUT_CONES=[[0,1,2,3]]);

julia> F=addStackStructure(X,[1,2,2,5]);

julia> showStackyFan(F);
```
"""   
function showStackyFan(F::StackyFan;stackypoints::Bool=true)
    X=F.fan
    plot()
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    dim=size(rayMatrix,2)
    if !(dim in [2,3])
         println("Fan visualization is only supported in dimensions 2 and 3.")
    end
    maxCones=convertIncidenceMatrix(X.MAXIMAL_CONES)
    for cone in maxCones
        if dim==2
            coneRays=slicematrix(rowMinors(rayMatrix,cone))
            plot2dCone(coneRays)
        elseif dim==3
            C=coneConvert(cone,rayMatrix)
            plot3dCone(coneVectorOrder(C))
        end
    end
    if stackypoints==true
        stackyPoints=slicematrix(findStackyRayMatrix(F))
        for point in stackyPoints
            if dim==2
                (x,y)=point
                plot!([x],[y],shape=:circle,color=:red)
            elseif dim==3
                (x,y,z)=point
                plot!([x],[y],[z],shape=:circle,color=:red)
            end
        end
    end
    plot!(legend=false)
end
