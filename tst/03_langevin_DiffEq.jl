using DifferentialEquations
using StochasticDelayDiffEq

using Statistics

using Plots
using FFMPEG


function pot(x;p=[0,0,0])
    a, b, K = p
    return a/4*x^4 - (b+K)/2*x^2 + (b+K)^2/(4*a)
end

function rhs_det(du, u, h, p, t)
    a, b, K, D0, τ = p
    du .= -a .* u.^3 .+ (b+K) .* u .- K .* h(p, t - τ)
end

function rhs_stoch(du, u, h, p, t)
    a, b, K, D0, τ = p
    du .= sqrt(2*D0)
end

h(p, t) = (ones(1) .+ t);

tmin =  0.0
tmax = 10.0
tpts = 301

tlen = tmax-tmin
dt = tlen/(tpts-1)

tspan = (0.0, 10.0)

par = [1.0, 10.0, 5.0, 0.1, 0.1]

rpts = 1000

Uvec = Vector{Vector{Float64}}(undef,rpts)
Tvec1 = Vector{Float64}(undef,rpts)

for ir in 1:rpts
    u0 = -0.01

    prob = SDDEProblem(
        rhs_det,
        rhs_stoch,
        [u0],
        h,
        (tmin,tmax),
        par;
        constant_lags = (par[end],)
    );

    sol = solve(
        prob,
        RKMil(),
        saveat = dt
    )

    Uvec[ir] = vec(reduce(hcat, sol.u))
    Tvec1 = sol.t
end

#plot(Tvec1,Uvec,legend = false)

xmin = -4.0
xmax =  4.0
xpts = 300

Xvec = collect(LinRange(xmin,xmax,xpts))
dx = (xmax-xmin)/(xpts-1)

rho = zeros((tpts,xpts)) 

for it in 1:tpts
    for ir in 1:rpts
        ix = Int(round((Uvec[ir][it]-xmin)/dx))+1
        rho[it,ix]+= 1.0/rpts
    end
end

rho = [rho[i,:] for i in 1:size(rho,1)]


video_path = "tst/output_video.mp4"
fps = 24

plot_ymax = mean(maximum.(rho)[2:end])


Vvec = pot.(Xvec,p=par)
Vvec.*=plot_ymax/maximum(Vvec)

anim = @animate for it in 1:tpts
    plot(Xvec,Vvec)
    plot!(Xvec,rho[it,:],legend=false,ylimits=(0,plot_ymax))
end

mp4(anim, video_path, fps=fps)

@info "done"