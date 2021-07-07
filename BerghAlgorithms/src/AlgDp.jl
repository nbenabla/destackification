export divAlongDivisor

"""
    divAlongDivisor(::Array{Int64,1},::StackyFan,::Dict)

    Calculates the simplified divisorial index along a divisor of a given cone in a fan with divisorial rays. Takes the subcone consisting of all relevant residual rays along with the specified divisor, and counts the number of residual rays that are relevant with respect to that subcone.

# Examples
```jldoctest AlgDp
julia> F=makeStackyFan([1 0 0; 0 1 0; 1 0 2],[[0,1,2]],[1,2,1])

julia> div=Dict([1,0,0]=>0,[0,1,0]=>0,[1,0,2]=>1)
    
julia> divAlongDivisor([1,2,3],F,div,[1,0,2])
2
julia> div=Dict([1,0,0]=>0,[0,1,0]=>1,[1,0,2]=>1) #[1 0 0] [1 0 2] - [1 0 0] has non-zero projection
    
julia> divAlongDivisor([1,2,3],F,div,[1,0,2])
1
julia> div=Dict([1,0,0]=>1,[0,1,0]=>0,[1,0,2]=>1) #[0 1 0] has stacky value >1
    
julia> divAlongDivisor([1,2,3],F,div,[1,0,2])
1
julia> div=Dict([1,0,0]=>0,[0,1,0]=>1,[1,0,2]=>1) #[1 0 0] [0 1 0] - [1 0 0] is independent
    
julia> divAlongDivisor([1,2,3],F,div,[0,1,0])
0
```
"""
function divAlongDivisor(cone::Array{Int64,1}, F::StackyFan, div::Dict, divisor::Array{Int64,1})
    # This may not be a very clear name, doesn't mention "simplified" or "index"
    X=deepcopy(F)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
    slicedRayMatrix=slicematrix(rayMatrix)
    
    # Construct the cone made up of the residual (non-divisorial) rays and the divisor
    residRays = Array{Int64, 1}[]
    residStack = Int64[]
    c = 0
    # Add the divisor (which should be divisorial)
    if div[divisor] == 1
        c += 1
        push!(residRays, divisor)
        push!(residStack, F.stacks[encode(divisor)])
    else
        error("The given ray is non-divisorial")
    end     
    # Append the other rays
    for i in cone
        ray = slicedRayMatrix[i]
        stack = F.stacks[encode(ray)]
        if div[ray] == 0 && isRelevant(ray, cone, F)
            c += 1
            push!(residRays, ray)
            push!(residStack, stack)
        end
    end
    residConeZ = Array(0:c-1)
    residCone = Array(1:c)
    # Create the stacky fan structure containing only the cone
    residRayMatrix = Array{Int64}(transpose(hcat(residRays...)))
    subconeFan = makeStackyFan(residRayMatrix, [residConeZ], residStack)
    
    # Count the number of residual rays which are dependent in the above cone
    # Skip the first ray, which is the divisor
    divIndex = 0
    for i in 2:c
        ray = residRays[i]
        if isRelevant(ray, residCone, subconeFan)
            divIndex += 1
        end
    end

    return divIndex
end

"""
    positiveDivIndexCone()

    Tests whether there exists a cone with non-zero divisorial index along the given divisor.
"""
function positiveDivIndexCone(divisor::Array{Int64,1}, F::StackyFan, div::Dict)
    # Go through each maximal cone containing the ray, and testing the divisorial index
    rayMatrix = convert(Array{Int64,2}, Array(Polymake.common.primitive(X.fan.RAYS)))
    maxconeList = convertIncidenceMatrix(F.fan.MAXIMAL_CONES)
    divisorIndex = getIndex(divisor, rayMatrix)
    for maxCone in maxconeList
        if divisorIndex in maxCone
            coneDivIndex = divAlongDivisor(maxCone, F, div, divisor)
            # If the divisorial index along the divisor is non-zero, return true
            if coneDivIndex > 0
                return true
            end
        end
    end
    return false
end
