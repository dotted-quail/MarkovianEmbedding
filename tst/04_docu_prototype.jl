using VoronoiFVM
using ExtendableGrids

using Plots
using FFMPEG

# Bernoulli function used in the exponential fitting discretization
function bernoulli(x)
    if abs(x) < nextfloat(eps(typeof(x)))
        return 1
    end
    return x / (exp(x) - 1)
end

function gauss(x,x0,s)
    g = 1/(sqrt(2*pi)*s)
    g*= exp(-1/2*((x-x0)/s)^2)
    return g
end


verbose=false
unknown_storage=:dense


struct data
    v::Vector
    D::Float64
    dx::Float64
end


xmin = -3.0
xmax =  3.0
xpts = 110
X = collect(range(xmin,stop=xmax,length=xpts))
dx = X[2]-X[1]

dat = data([1.0],0.1,dx)

grid = VoronoiFVM.Grid(X)

function force(x,data)
    dx = data.dx
    return dat.v * 2 * dx * (1.0 - dx), 0
end

dum = Vector{Any}(undef,20)

function exponential_flux!(f,u,edge,data)
    vh = project(edge, data.v)
    #global dum[1] = edge
    #global dum[2] = data.v
    #global dum[3] = vh
    Bplus = data.D * bernoulli(vh / data.D)
    Bminus = data.D * bernoulli(-vh / data.D)
    f[1] = Bminus * u[1,1] - Bplus * u[1,2]
end


function storage!(f,u,node,data)
    f[1] = u[1]
end

physics = VoronoiFVM.Physics(
    flux = exponential_flux!,
    storage = storage!,
    data = dat
)

sys = VoronoiFVM.System(
    grid,
    physics,
    unknown_storage = unknown_storage
)

# Add species 1 to region 1
enable_species!(sys,1,[1])


# Create a solution array
inival=unknowns(sys)
solution=unknowns(sys)

# Broadcast the initial value
x0 = 0.0
s = 0.1
inival[1,:].=map(x->gauss(x,x0,s),X)


# Create solver control info
control=VoronoiFVM.NewtonControl()
control.verbose=verbose


tmin = 0.0
tmax = 5
tpts = 400
T = collect(range(tmin,stop=tmax,length=tpts))
dt = T[2]-T[1]

sol_arr = Array{Float64}(undef, tpts, xpts)



for it in 1:tpts
    solve!(solution,inival,sys,control=control,tstep=dt)
    inival .= solution
    sol_arr[it,:] = vec(solution)
end


video_path = "tst/output_video.gif"
fps = 24

#=
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
end

gif(anim, video_path, fps=fps)
=#
