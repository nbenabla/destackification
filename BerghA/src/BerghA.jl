module BerghA
include("StackyFan.jl")
export BerghAEfficient

# Write your package code here.
"""

    BerghAEfficient(F::StackyFan,D::Array{Int64,1})

    Given a stacky fan F and a vector of booleans D representing the distinguished structure,
    returns a smooth stacky fan where the distinguished rays are independent.

    
"""

function BerghAEfficient(F::StackyFan,D::Array{Int64,1};verbose::Bool=false)
    if verbose==true
        println("==algorithm is running in verbose mode==")
        println(" ")
        println("=======")
    end
    X=deepcopy(F)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
    coneList=getCones(X.fan)
    dim=size(rayMatrix,2)
    numRays=size(rayMatrix,1)
    
    #check if the vector D has length equal to the number of rays in F
    if numRays != size(D,1)
        error("length of vector representing distinguished structure does not agree with number of rays in stacky fan.")
    end
    
    #A0: initialization
    i=0
    while(true)
        rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
        numRays=size(rayMatrix,1)
        coneList=getCones(X.fan)
        
       #Debugging
        #if verbose==true
            #println("Ray matrix: $rayMatrix")
            #println("Cone list: $coneList")
            #sW=stackyWeights(X)
            #println("Stacky weights: $sW")
            #println("Distinguished rays: $D")
        #end
        #coneMultiplicities=Int64[]
        #for cone in coneList
        #    C=coneConvert(cone,rayMatrix)
        #    push!(coneMultiplicities,coneMultiplicity(C))
        #end
        #if verbose==true
            #println("Cone multiplicities: $coneMultiplicities")
        #end
       #End debugging
        
        #A1: Find S the set of cones that contain a distinguised ray and an interior lattice point 
        #Note: cones in S are 1-indexed.
        S=filter(cone->distinguishedAndIntPoint(cone,rayMatrix,D),coneList)
        # If S is empty, the program terminates.
        if S==[]
            break
        end
        
        #A2 - find extremal cones
        Smax=extremalCones(S,rayMatrix,D)
        
        #Print information on the number of extremal cones, their number of non-distinguished rays, and their multiplicity
        #The algorithm is structured to first reduce the number of non-distinguished rays in extremal cones, and then reduce the multiplicity of said cones,
            #so this information can be used to track the algorithm's progress
        if verbose==true
            Smaxcount=size(Smax,1)
            println("Number of extremal cones: $Smaxcount")
            testCone=Smax[1]
            c1=convertToIncidence(testCone,numRays)
            nonDist=size(testCone,1)-dot(c1,D)
            mult=coneMultiplicity(coneConvert(testCone,rayMatrix))
            println("Maximal non-distinguished rays and multiplicity: $nonDist, $mult")
        end
        
        #A2 - find interior points in Smax
        intPoints=[]
        for cone in Smax
            #C=rowMinors(rayMatrix,cone)
            C=coneConvert(cone,rayMatrix)
            coneIntPoints=interiorPoints(C)
            for point in coneIntPoints
               push!(intPoints,(point,cone)) #the point is stored as a tuple along with the cone containing it
            end
        end
        
        #A2 - find stacky points (in terms of coefficients) derived from interior points
        P=Array{Int64,1}[]
        for (point,cone) in intPoints
            stackyPoint=coneRayDecomposition(cone,rayMatrix,point,stackyWeights(X)) #each interior point is rewritten as a string of coefficients
                #corresponding to its representation as a sum of stacky rays
            push!(P,stackyPoint)
        end
        
        #A2 - find element of P such that the sum of the entries is minimal.
            #This invariant does not produce a unique ray, so there is a degree of arbitrary selection.
        psi=minimalByDist(P,D)
        #if verbose==true
            #println("Psi: $psi")
        #end
        
        
        #A3 - perform root construction
        X=rootConstructionDistinguishedIndices(X,D,psi)
        
        #A3 - modify psi with respect to root construction
        for i in 1:length(psi)
            if D[i]==1 && psi[i]>0
                psi[i]=1
            end
        end
        
        #A5 - perform repeated stacky barycentric star subdivision with respect to psi.
        while(count(x->x>0,psi)>1)
            #A4 - perform stacky star subdivision
            # Get the indices of the non-zero coefficients in psi - this is used to determine the cone 
                #containing the support of psi, which will be subdivided
            
            supportCone=findall(x->x!=0,psi)
            #find the stacky barycenter of that cone, which becomes the exceptional (blowup) ray
            exceptional=findStackyBarycenter(supportCone,X)
            
            #Performing a blowup may cause the rays defined in the StackyFan struct to be reordered - this property is inherited from
                #the Oscar/Polymake fan object. Since D (distinguished rays) and psi are lists, they must be reordered to match the
                #new order of the rays of X. This reordering is accomplished through defining dictionaries before the blowup is performed.
            
            code_rays = mapslices(encode, Polymake.common.primitive(X.fan.RAYS), dims=2)
            # Track the indices of distinguished rays
            D_pairs = map((x,y) -> (x,y), code_rays, D)
            D_Dict = Dict(D_pairs)
            # Track psi as a linear combination of the generators
            psiPairs = map((x,y) -> (x,y), code_rays,psi)
            psiDict = Dict(psiPairs)

            #perform the blowup
            X=stackyBlowup(X,[x-1 for x in supportCone],exceptional)
            
            G=gcd(exceptional) #since the blowup ray may not be primitive, it is made primitive and then assigned a stacky value so its stacky form is unchanged.
            primExcep=Polymake.common.primitive(exceptional)
            
            # Update the dictionaries storing fan information
            D_Dict[encode(primExcep)]=1
            psiDict[encode(primExcep)]=1
            
            #create new lists
            newRays=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS))))
            newD=Int64[]
            newpsi=Int64[]
            for ray in newRays
                E=encode(ray)
                excepCode=encode(primExcep)
                push!(newD,D_Dict[E])
                #A4 - modify psi
                if E==excepCode
                    push!(newpsi,1)
                elseif psiDict[E]>1
                    push!(newpsi,psiDict[E]-1)
                else
                    push!(newpsi,0)
                end
            end
            psi=newpsi
            D=newD
        end
        if verbose==true
            println("=======")
        end
        i+=1
    end
    if verbose==true
        println("Number of steps: $i")
    end
    return X
end

end
