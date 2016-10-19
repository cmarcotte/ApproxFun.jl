using ApproxFun, Base.Test
import Compat: view





## Try constructor variants

ff=(x,y)->exp(-10(x+.2)^2-20(y-.1)^2)*cos(x*y)
gg=x->exp(-10(x[1]+.2)^2-20(x[1]-.1)^2)*cos(x[1]*x[2])
f=Fun(ff,Interval()^2,10000)
@test_approx_eq f(0.,0.) ff(0.,0.)

f=Fun(gg,Interval()^2,10000)
@test_approx_eq f(0.,0.) ff(0.,0.)

f=Fun(ff,Interval()^2)
@test_approx_eq f(0.,0.) ff(0.,0.)
f=Fun(gg,Interval()^2)
@test_approx_eq f(0.,0.) ff(0.,0.)


f=Fun(ff)
@test_approx_eq f(0.,0.) ff(0.,0.)
f=Fun(gg)
@test_approx_eq f(0.,0.) ff(0.,0.)


# Fun +-* constant
f=Fun((x,y)->exp(x)*cos(y))

@test_approx_eq f(0.1,0.2)+2 (f+2)(0.1,0.2)
@test_approx_eq f(0.1,0.2)-2 (f-2)(0.1,0.2)
@test_approx_eq f(0.1,0.2)*2 (f*2)(0.1,0.2)


## LowRankFun

F = LowRankFun((x,y)->besselj0(10(y-x)),Chebyshev(),Chebyshev())

@test_approx_eq F(.123,.456) besselj0(10(.456-.123))

G = LowRankFun((x,y)->besselj0(10(y-x));method=:Cholesky)

@test_approx_eq G(.357,.246) besselj0(10(.246-.357))





## 1D in 2D

d=Interval((0.,0.),(1.,1.))
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))
@test_approx_eq f(FixedSizeArrays.Vec(0.5,0.5)) exp(-0.5-2cos(0.5))

f=Fun(xy->exp(-xy[1]-2cos(xy[2])),d,20)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))

f=Fun((x,y)->exp(-x-2cos(y)),d,20)
@test_approx_eq f(0.5,0.5) exp(-0.5-2cos(0.5))


d=Circle((0.,0.),1.)
f=Fun(xy->exp(-xy[1]-2cos(xy[2])),Fourier(d),40)
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))
@test_approx_eq f(FixedSizeArrays.Vec(cos(0.1),sin(0.1))) exp(-cos(0.1)-2cos(sin(0.1)))

f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d),40)
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))


f=Fun((x,y)->exp(-x-2cos(y)),Fourier(d))
@test_approx_eq f(cos(0.1),sin(0.1)) exp(-cos(0.1)-2cos(sin(0.1)))

println("    Calculus tests")

## Sum

ff=(x,y)->(x-y)^2*exp(-x^2/2.-y^2/2)
f=Fun(ff,Domain([-4.,4.])^2)
@test_approx_eq f(0.1,0.2) ff(0.1,0.2)

@test_approx_eq sum(f,1)(0.1) 2.5162377980828357
f=LowRankFun(f)
@test_approx_eq evaluate(f.A,0.1) map(f->f(0.1),f.A)


## Kron operator
Mx=Multiplication(Fun(cos),Chebyshev())
My=Multiplication(Fun(sin),Chebyshev())
K=Mx⊗My
@test_approx_eq ApproxFun.BandedBlockBandedMatrix(view(K,1:10,1:10)) K[1:10,1:10]
C=Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical(1)⊗Ultraspherical(1))
@test_approx_eq C[1:100,1:100] Float64[C[k,j] for k=1:100,j=1:100]


# 2d derivative (issue #346)
let d = Chebyshev()^2
    f = Fun((x,y) -> sin(x) * cos(y), d)
    C=Conversion(Chebyshev()⊗Chebyshev(),Ultraspherical(1)⊗Ultraspherical(1))
    @test_approx_eq (C*f)(0.1,0.2) f(0.1,0.2)
    Dx = Derivative(d, [1,0])
    f = Fun((x,y) -> sin(x) * cos(y), d)
    fx = Fun((x,y) -> cos(x) * cos(y), d)
    @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
    Dy = Derivative(d, [0,1])
    fy = Fun((x,y) -> -sin(x) * sin(y), d)
    @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
    L=Dx+Dy
    @test_approx_eq (L*f)(0.2,0.3) (fx(0.2,0.3)+fy(0.2,0.3))

    B=ldirichlet(d[1])⊗ldirichlet(d[2])
    @test_approx_eq B*f f(-1.,-1.)

    B=Evaluation(d[1],0.1)⊗ldirichlet(d[2])
    @test_approx_eq B*f f(0.1,-1.)

    B=Evaluation(d[1],0.1)⊗Evaluation(d[2],0.3)
    @test_approx_eq B*f f(0.1,0.3)

    B=Evaluation(d,(0.1,0.3))
    @test_approx_eq B*f f(0.1,0.3)
end

let d = Space([0,1]) * Space([0,2])
    Dx = Derivative(d, [1,0])
    f = Fun((x,y) -> sin(x) * cos(y), d)
    fx = Fun((x,y) -> cos(x) * cos(y), d)
    @test (Dx*f)(0.2,0.3) ≈ fx(0.2,0.3)
    Dy = Derivative(d, [0,1])
    fy = Fun((x,y) -> -sin(x) * sin(y), d)
    @test (Dy*f)(0.2,0.3) ≈ fy(0.2,0.3)
    L=Dx+Dy
    @test_approx_eq (L*f)(0.2,0.3) (fx(0.2,0.3)+fy(0.2,0.3))

    B=ldirichlet(d[1])⊗ldirichlet(d[2])
    @test abs(B*f-f(0.,0.)) ≤ 10eps()

    B=Evaluation(d[1],0.1)⊗ldirichlet(d[2])
    @test_approx_eq B*f f(0.1,0.)

    B=Evaluation(d[1],0.1)⊗Evaluation(d[2],0.3)
    @test_approx_eq B*f f(0.1,0.3)

    B=Evaluation(d,(0.1,0.3))
    @test_approx_eq B*f f(0.1,0.3)
end



## x,y constructor


x,y=Fun(Interval()^2)
@test_approx_eq x(0.1,0.2) 0.1
@test_approx_eq y(0.1,0.2) 0.2
