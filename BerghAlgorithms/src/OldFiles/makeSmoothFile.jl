using Oscar
using Combinatorics
export makeSimplicial
export makeSmooth

function printn(A)
    print(A)
    print("\n")
end


"""
    slicematrix(::AbstractMatrix{<:Number})

    Take a two-dimensional matrix and output a list of its row vectors.

# Examples
```jldoctest makeSmoothWithDependencies
julia> A=[1 2; 3 4]

julia> slicematrix(A)
[[ 1 ,  2 ], [ 3 ,  4 ]]
"""

function slicematrix(A::AbstractMatrix{<:Number})
    return [A[i, :] for i in 1:size(A,1)]
end

"""
    rowMinors(::AbstractMatrix{<:Number},::Union{AbstractSet,AbstractVector})

    Identical to slicematrix, except only returns row vectors indexed by a set S.

# Examples
```jldoctest makeSmoothWithDependencies
julia> A=[1 2 3;4 5 6; 7 8 9]

julia> S=Set([1,3])

julia> rowMinors(A,S)
2×3 LinearAlgebra.Transpose{Int64,Array{Int64,2}}:
 1  2  3
 7  8  9
"""

function rowMinors(A::AbstractMatrix{<:Number},S::Union{AbstractSet,AbstractVector})
    outList=[]
    slices=slicematrix(A)
    for i in 1:size(slices,1)
        if i in S
            append!(outList,[slices[i]])
        end
    end
    return transpose(hcat(outList...))
end

"""
    convertIncidenceMatrix(::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})

    Takes a Polymake incidence matrix (e.g., the output of X.MAXIMAL_CONES for a toric variety X) and outputs a list of vectors,
    with each vector recording the indices marked on a given row of the incidence matrix.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]])

julia> M=X.MAXIMAL_CONES

julia> convertIncidenceMatrix(M)
[[ 1 ,  2 ], [ 2 ,  3 ]]
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
    return convert.(Array{Int64, 1}, out)
end

"""

    findBarycenter(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated)

    Takes a normal toric variety X and a set s corresponding to a subset of rays of X, and outputs a polymake vector
    corresponding to the barycenter of those rays.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]])

julia> s=[1,2]

julia> findBarycenter(s,X)
pm::Matrix<pm::Integer>
2 1

"""

function findBarycenter(s::Union{AbstractSet,AbstractVector},X::Polymake.BigObjectAllocated)
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

    coneConvert(::abstractVector{Int64},::abstractMatrix{Int64})

    Takes a matrix where the columns represent rays, and a list of indices, and forms a Polymake cone object generated by the rays corresponding to those indices.

# Examples
```jldoctest StackyFan

julia> typeof(coneConvert([1, 2, 4],[1 0 0; 0 1 0; 0 0 1; 1 1 1]))
Polymake.BigObjectAllocated

"""

function coneConvert(cone::Array{Int64,1},rayMatrix::Array)
    coneRays=rowMinors(rayMatrix,cone)
    C=Polymake.polytope.Cone(RAYS=coneRays)
    return C
end

"""

    getCones(X::Polymake.BigObjectAllocated)
    
    Returns all the cones of a fan X as a list of lists, with each interior list containing the indices of the rays generating a given cone.

# Examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 ],INPUT_CONES=[[0,1,2,3]])

julia> getCones(X)
[[ 0 ,  1 ,  2 ,  3 ], [ 0 ,  1 ], [ 0 ,  2 ], [ 2 ,  3 ], [ 1 ,  3 ], [ 0 ], [ 1 ], [ 2 ], [ 3 ]]



"""

function getCones(X::Polymake.BigObjectAllocated)
    lattice=X.HASSE_DIAGRAM
    faces=@Polymake.convert_to Array{Set{Int}} lattice.FACES
    out=[]
    for i in 2:(size(faces,1)-1)
        push!(out,Array(@Polymake.convert_to Array{Int} faces[i]))
    end
    return out
end

function findFaceContainingRay(C::Polymake.BigObjectAllocated,v::Array{Int64, 1})
    rayMatrix=Array(C.RAYS)
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

function findMinimalCone(X::Polymake.BigObjectAllocated,v::Array{Int64, 1})
    rayMatrix=Array(X.RAYS)
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
    

    
function starSubdivision(X::Polymake.BigObjectAllocated, v::Array{Int64, 1})
    minimalCone = findMinimalCone(X, v)
    s = [i - 1 for i in minimalCone]
    v = transpose(v)
    return toric_blowup(s, X, v)
end

"""

    toric_blowup(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated,::AbstractVector)

    Takes a normal toric variety X, a set s corresponding to a subset of rays of X, and a (optional) polymake vector,
    v, blow up X at v. If v is not provided, blow up X at the barycenter of s.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2,3]])

julia>toric_blowup([0, 1, 2, 3], X, nothing)
INPUT_CONES
{0 1 4}
{0 2 4}
{2 3 4}
{1 3 4}
INPUT_RAYS
1 0 0
1 1 0
1 0 1
1 1 1
2 1 1

julia>toric_blowup([0], X, nothing)
INPUT_CONES
{0 2 3}
{0 1 3}
INPUT_RAYS
1 0 0
1 1 0
1 0 1
1 1 1

"""

function toric_blowup(s, X, v)
    s = [i + 1 for i in s]
    # If v is not provided, take v as the barycenter
    if v==nothing
        v=findBarycenter(s,X)
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
    # Remove duplicates
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

    convertBool(::AbstractVector)

    Takes a column vector of boolean values and converts it to a vector of indices marked 'true'.

#Examples
```jldoctest makeSmoothWithDependencies
julia> B=[true true false true]

julia> convertBool(transpose(B))
[0, 1, 3]
"""

function convertBool(B::AbstractVector)
    out=[]
    for i in 1:size(B,1)
        if B[i]==true
           append!(out,i-1) 
        end
    end
    return out
end


"""

    getConeRank(::AbstractMatrix,::AbstractVector)

    Takes a matrix and a vector containing indices corresponding to rows of a matrix,
    and calculates the rank of the matrix consisting only of those rows.

#Examples
```jldoctest makeSmoothWithDependencies
julia> v=[1,2]

julia> M=[0 1; 1 1; 1 0]

julia> getConeRank(v,M)
2
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
    makeSimplicial(::Polymake.BigObjectAllocated)

    Takes in a normal toric variety and returns a simplicial toric variety 
    by subdividing (blowing up) the non-simplicial maximal cones.

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

    Takes in a normal toric variety X and output a new smooth toric variety by iteratively blowing up.
    In the language of fans, these blowups are achieved by subdividing non-smooth cones.

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

"""

function makeSmooth(X::Polymake.BigObjectAllocated)
    Y  = copy(X)
    while(true)
        coneList = convertIncidenceMatrix(Y.MAXIMAL_CONES)
        rayMatrix = Array(Y.RAYS)
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