using MarkovianEmbedding

using DifferentialEquations
using JSON
using Plots
using Random


function rhs_det(y, p)
    a = p.a
    b = p.b
    K = p.K

    return -a*y^3+(b+K)*y
end

function rhs_del(yt, p)
    K = p.K

    return -a*y^3+(b+K)*y
end

function rhs_stoch(y, p)
    D0 = p.D0
    rnd = randn(p.rng, Float64)

    return sqrt(2*D0)*rnd
end

function pot(x, p)
    a = p.a
    b = p.b

    return a/4*x^4 - b/2*x^2 + b^2/(4*a)
end

mutable struct par
    a::Float64
    b::Float64
    K::Float64

    D0::Float64
    rng::AbstractRNG
    
    tau::Float64
end

a  = 1.0
b  = 10.0
K  = 0.0
D0 = 0.001
rng = MersenneTwister(1234)
tau = 0.1

P = par(a,b,K,D0,rng,tau)

xmin = -5
xmax =  5
xpts =  2^9
Xvec = collect(LinRange(xmin,xmax,xpts))

V = P.a/4*Xvec.^4-P.b/2*Xvec.^2 .+ P.b^2/(4*P.a)

y0 = 3.0

tmin = 0.0
tmax = 5.0
tpts = 500
Tvec = collect(range(tmin,stop=tmax,length=tpts+1)[1:end-1])
dt = Tvec[2]-Tvec[1]

taupts = 4


Yvec = Vector{Float64}(undef,tpts)
Yvec[1] = y0

for it in 2:tpts
    y = Yvec[it-1]
    Yvec[it] = y
    Yvec[it]+= rhs_det(y, P)*dt
    Yvec[it]+= rhs_stoch(y, P)*sqrt(dt)
end


#plot(Xvec,V,label='V')
#scatter!([x0],[pot(x0,P)])

plot(Tvec,Yvec)






