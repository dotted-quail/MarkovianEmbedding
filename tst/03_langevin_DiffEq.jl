using DifferentialEquations
using StochasticDelayDiffEq

using Statistics

using Plots
using FFMPEG

using JSON


par = [1.0, 10.0, 1.0, 0.1, 0.1]

begin
    a   = 1.0
    b   = 10.0
    K   = 1.0
    D0  = 0.1
    tau = 0.1
end

json_path = pwd()*"/tst/"
json_name = "test_par.json"
json_dir = json_path*json_name

system_parameters = Dict(
    "a" => a,
    "b" => b,
    "K" => K,
    "D0" => D0,
    "tau" => tau
)

output_dict = Dict(
    "system_parameters" => system_parameters
)

open(json_dir, "w") do file
    JSON.print(file, output_dict)
end

open(json_dir, "r") do file
    input_dict = JSON.parse(file)
end













function pot(x;p=[0,0,0,0,0])
    a, b, K, D0, τ = p
    K=0
    return a/4*x^4 - (b)/2*x^2 + (b)^2/(4*a)
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

xmin = -6.0
xmax =  6.0
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