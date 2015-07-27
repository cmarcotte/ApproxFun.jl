# A DiracSpace represents dirac delta functions.
# A PointSpace represents finite functions defined on points


export deltafunction

abstract AbstractPointSpace <:FiniteSpace{RealBasis} #TODO: SingularBasis


for TYP in (:DiracSpace,:PointSpace)
    @eval begin
        immutable $TYP{T<:Number}<:AbstractPointSpace
            points::Vector{T}
        end

#         $TYP{T<:Number}(d::Vector{T})=$TYP{T}(d)
#         $TYP{T<:Integer}(d::Vector{T})=$TYP{Float64}(d)
#         $TYP{P<:Point}(d::Vector{P})=$TYP(map(p->p.x,d))

        $TYP() = $TYP(Float64[])
        setdomain(DS::$TYP,d::UnionDomain)=$TYP(d)

        spacescompatible(a::$TYP,b::$TYP)= a.points==b.points

        union_rule(a::$TYP,b::$TYP)=$TYP(sort(union(a.points,b.points)))
    end
end

Base.length(C::AbstractPointSpace)=length(C.points)


domain(DS::AbstractPointSpace)=mapreduce(Point,union,DS.points)

canonicalspace(a::AbstractPointSpace)=a




# DiracSpace specific
function Base.sum{DS<:DiracSpace}(f::Fun{DS})
    n = length(space(f).points)
    sum(f.coefficients[1:n])
end


function evaluate{DS<:DiracSpace}(f::Fun{DS},x::Number)
  n = length(f.space.points)
  if x in f.space.points
    error("You cannot evaluate a Dirac delta at its center.")
  else
    evaluate(Fun(f.coefficients[n+1:end],f.space.space),x)
  end
end

evaluate(S::DiracSpace,coeffs::Vector,x::Number)=evaluate(Fun(coeffs,S),x)



function coefficients(cfs::Vector,fromspace::DiracSpace,tospace::DiracSpace)
    if spacescompatible(fromspace,tospace)
        return cfs
    end


    cfs=pad(cfs,max(length(tospace.points),length(fromspace.points)))

    # this first for-loop removes coefficients of Dirac points that are zero
    nonzerofromspacepoints = Float64[]
    nonzeroDiraccfs = Float64[]
    for i = 1:length(fromspace.points)
        if cfs[i] != 0
            push!(nonzerofromspacepoints, fromspace.points[i])
            push!(nonzeroDiraccfs, cfs[i])
        end
    end

    # if the points that remain can be represented in the tospace
    if issubset(nonzerofromspacepoints,tospace.points)
        finalDiraccfs = []
        j=1 #counter for the nonzerofromspacepoints
        for i = 1:length(tospace.points)
            if nonzerofromspacepoints[j] in tospace.points
                finalDiraccfs[end] = nonzeroDirccfs[j]
                j += 1
            else
                finalDiraccfs[end]=0
            end
        end
        finalDiraccfs
    else
        error("The space you are converting from has Dirac deltas that cannot be represented in the space you are converting to.")
    end
end




## deltafunction

deltafunction(x)=Fun([1.0],DiracSpace([x]))

