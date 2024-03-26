using VoronoiFVM
using ExtendableGrids

using Plots
using FFMPEG

struct data
    a::Float64
    b::Float64
    K::Float64
    D0::Float64
    s::Float64
    x0::Float64
end

dat = data(1,10,1,0.1,1.0,0.0)

xmin = -6.0
xmax =  6.0
xpts = 110
X = collect(range(xmin,stop=xmax,length=xpts))
dx = X[2]-X[1]

grid = simplexgrid(X)

tmin = 0.0
tmax = 0.7
tpts = 400
T = collect(range(tmin,stop=tmax,length=tpts))
dt = T[2]-T[1]



function force(x,n)
    a = 1.0
    b = 10.0
    K = 1.0
    return (-(a*x^3-x*b-x*K),0)
end

evelo = edgevelocities(grid, force)

function flux!(f,u,edge,data)
    D0 = data.D0
    vd = evelo[edge.index] / D0
    bp = fbernoulli(vd)
    bm = fbernoulli(-vd)
    f[1] = D0 * (bp * u[1] - bm * u[2])
end

function storage!(f,u,node,data)
    f[1] = u[1]
end

phy = VoronoiFVM.Physics(
    flux = flux!,
    storage = storage!,
    data = dat
)
sys = VoronoiFVM.System(grid,phy)

enable_species!(sys, 1, [1])

function gauss(x,data)
    s = data.s
    x0 = data.x0
    g = 1/(sqrt(2*pi)*s)
    g*= exp(-1/2*((x-x0)/s)^2)
    return g
end



inival = unknowns(sys)
sol =  unknowns(sys)
inival[1,:] .= map(x->gauss(x,dat),X)

sol_arr = Array{Float64}(undef, tpts, xpts)

for it in 1:tpts
    solve!(sol,inival,sys,tstep=dt)
    inival .= sol
    sol_arr[it,:] = vec(sol)
end

video_path = "tst/output_video.gif"
fps = 24


V = dat.a/4*X.^4 - (dat.b+dat.K)/2*X.^2 .+ (dat.b+dat.K)^2/(4*dat.a)


anim = @animate for it in 1:tpts
    Plots.plot(
        X,sol_arr[it,:],
        label = "",
        xlimits=(xmin,xmax),
        ylimits=(minimum(sol_arr),maximum(sol_arr)),
        xlabel = "\$x\$",
        ylabel = "\$\\rho_1\$",
        size = (800,400),
        margin=5Plots.mm,
        title="t="*string(round(T[it],sigdigits=3))
    )
    Plots.plot!(X,V*0.1)
end

gif(anim, video_path, fps=fps)