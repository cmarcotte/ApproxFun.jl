using ApproxFun, Base.Test, Compat, Base.Test
    import ApproxFun: bandedblockbandedoperatortest, bandedblockoperatortest

## Check operators
S=JacobiWeight(1.,1.,Jacobi(1.,1.))^2
Δ=Laplacian(S)

bandedblockbandedoperatortest(Δ)

u=Fun((x,y)->sin(π*x)*sin(π*y),S)
f=-2π^2*u


QR=qrfact(Δ)
v=QR\f
@test norm((u-v).coefficients)<100eps()

v=Δ\f
@test norm((u-v).coefficients)<100eps()


f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),rangespace(Δ))  #default is [-1,1]^2
@time v=linsolve(Δ,f;tolerance=1E-14)
@test norm((Δ*v-f).coefficients)<1E-14



## Rectangle PDE
dx=dy=Interval()
d=dx*dy
g=Fun((x,y)->exp(x)*cos(y),∂(d))

bandedblockbandedoperatortest(Laplacian(d))
bandedblockoperatortest(Dirichlet(d))

A=[Dirichlet(d);Laplacian(d)]

u=A\[g,0.]
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))

bandedblockbandedoperatortest(Laplacian(d)+0.0I)

A=[Dirichlet(d);Laplacian(d)+0.0I]
u=A\[g,0.]
@test_approx_eq u(.1,.2) real(exp(0.1+0.2im))


println("    Poisson tests")


## Poisson

f=Fun((x,y)->exp(-10(x+.2)^2-20(y-.1)^2),Interval()^2,500)  #default is [-1,1]^2
d=domain(f)
A=[Dirichlet(d);Laplacian(d)]
@time  u=linsolve(A,[zeros(∂(d));f];tolerance=1E-7)
@test_approx_eq_eps u(.1,.2) -0.04251891975068446 1E-5


println("    Periodic Poisson tests")



d=PeriodicInterval()^2
f=Fun((x,y)->exp(-10(sin(x/2)^2+sin(y/2)^2)),d)
A=Laplacian(d)+.1I
u=A\f
@test (lap(u)+.1u-f)|>coefficients|>norm < 1000000eps()





# fourth order

println("    Bilaplacian tests")
dx=dy=Interval()
d=dx*dy
Dx=Derivative(dx);Dy=Derivative(dy)
L=Dx^4⊗I+2*Dx^2⊗Dy^2+I⊗Dy^4

bandedblockbandedoperatortest(L)

A=[dirichlet(dx)⊗eye(dy);
        eye(dx)⊗dirichlet(dy);
        neumann(dx)⊗eye(dy);
        eye(dx)⊗neumann(dy);
         L]


u=linsolve(A,ones(4);tolerance=1E-5)
@test_approx_eq u(0.1,0.2) 1.0


F=[Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[1]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[2]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[3]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[4]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[5]));
    Fun((x,y)->real(exp(x+1.0im*y)),rangespace(A[6]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[7]));
    Fun((x,y)->-imag(exp(x+1.0im*y)),rangespace(A[8]));
    0]

u=linsolve(A,F;tolerance=1E-10)

@test_approx_eq u(0.1,0.2)  exp(0.1)*cos(0.2)



## Test periodic x interval

println("    Periodic x Interval tests")
d=PeriodicInterval()*Interval()

u_ex=Fun((x,y)->real(cos(x+im*y)),d)

B=Dirichlet(Space(d))

g=Fun((x,y)->real(cos(x+im*y)),rangespace(B))  # boundary data
@test norm((B*u_ex-g).coefficients) < 10eps()

bandedblockbandedoperatortest(Laplacian(d))

u=[B;Laplacian(d)]\[g;0.]

@test_approx_eq u(.1,.2) real(cos(.1+.2im))


dθ=PeriodicInterval(-2.,2.);dt=Interval(0,1.)
d=dθ*dt
Dθ=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
B=eye(d[1])⊗ldirichlet(dt)

u0=Fun(θ->exp(-20θ^2),dθ,20)

bandedblockbandedoperatortest(Dt+Dθ)
u=linsolve([B;Dt+Dθ],[u0;0.];tolerance=1E-7)

@test_approx_eq_eps u(.1,.2) u0(0.1-0.2) 1E-7



## Small diffusoion


println("    Time evolution tests")


dx=Interval();dt=Interval(0,0.2)
d=dx*dt
Dx=Derivative(d,[1,0]);Dt=Derivative(d,[0,1])
x,t=Fun(dx*dt)


B=0.0
C=0.0
V=B+C*x
ε=0.1
f=Fun(x->exp(-30x^2),dx)

bandedblockbandedoperatortest(Dt-ε*Dx^2-V*Dx)

u=linsolve([timedirichlet(d);Dt-ε*Dx^2-V*Dx],[f;zeros(3)];tolerance=1E-6)

@test_approx_eq u(.1,.2) 0.496524222625512


## Schrodinger

dx=Interval(0.,1.);dt=Interval(0.0,0.001)
C=Conversion(Chebyshev(dx)*Ultraspherical(1,dt),Ultraspherical(2,dx)*Ultraspherical(1,dt))
bandedblockbandedoperatortest(C)
bandedblockbandedoperatortest(Operator{Complex128}(C))


d=dx*dt

x,y=Fun(d)
V=x^2

Dt=Derivative(d,[0,1]);Dx=Derivative(d,[1,0])

ϵ=1.
u0=Fun(x->exp(-100*(x-.5)^2)*exp(-1./(5*ϵ)*log(2cosh(5*(x-.5)))),dx)
L=ϵ*Dt+(.5im*ϵ^2*Dx^2)
bandedblockbandedoperatortest(L)

@time u=linsolve([timedirichlet(d);L],[u0;zeros(3)];tolerance=1E-5)
@test_approx_eq u(0.5,0.001) 0.857215539785593+0.08694948835021317im  # empircal from schurfact


## Periodic

println("    Periodic tests")

d=PeriodicInterval()^2
f=Fun((θ,ϕ)->exp(-10(sin(θ/2)^2+sin(ϕ/2)^2)),d)
A=Laplacian(d)+.1I
bandedblockbandedoperatortest(A)

@time u=A\f
@test_approx_eq u(.1,.2) u(.2,.1)
