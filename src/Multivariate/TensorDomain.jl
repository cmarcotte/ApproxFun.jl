

type TensorDomain{D<:Domain}
    domains::Vector{D} 
end

TensorDomain(A,B)=TensorDomain([A,B])
⊗(A::Domain,B::Domain)=TensorDomain(A,B)

domain(f::Fun2D)=domain(f.A[1])⊗domain(f.B[1])
