using MarkovianEmbedding

#=
a  = 1.0
b  = 0.5
K  = 0.1
D0 = 0.0001

tau = 1.0

tmin =  0.0
tmax = 10.0
tpts = 500

tlen = tmax-tmin
dt   = tlen/tpts
taupts = Int(round(tau/dt))

xmin = -5.0
xmax =  5.0
xpts =  2^7
Xvec = collect(LinRange(xmin,xmax,xpts))

xlen = xmax-xmin
dx   = xlen/xpts

ymin = -5.0
ymax =  5.0
ypts =  2^7
Yvec = collect(LinRange(ymin,ymax,ypts))

ylen = ymax-ymin
dy   = ylen/ypts

s0 = 0.01
x0 = 0.0

p1 = 1/sqrt(2*pi*s0^2)*exp.(-(Xvec.-x0).^2/(2*s0^2))
p2 = zeros((xpts,ypts))

p1_hist = Vector{Vector{Float64}}(undef,taupts)
p2_hist = Vector{Matrix{Float64}}(undef,taupts)

for itau in 1:taupts
    p1_hist[itau] = p1
    p2_hist[itau] = p2
end

write_cond = true
res_p1 = Vector{Vector{Float64}}(undef,tpts)
#res_p2 = 

for it in 1:tpts-1
    p1new = Vector{Float64}(undef,xpts)
    for ix in 2:xpts-1
        drif = ((3*a*Xvec[ix+1]^2-b-K)*p1[ix+1] - (3*a*Xvec[ix-1]^2-b-K) * p1[ix-1]) / (2*dx)
        diff = D0*(p1[ix+1] - 2*p1[ix] + p1[ix-1]) / dx^2

        p1new[ix] = p1[ix] - dt * (drif + diff)
    end

    #neumann conditions
    p1new[1] = p1new[2]
    p1new[end] = p1new[end-1]

    global p1 = p1new

    if write_cond
        res_p1[it] = p1new
    end
end



using Plots

# Example data (replace with your actual data)
vec_of_vecs = res_p1
vec_of_mats = [rand(10, 10) for _ in 1:length(vec_of_vecs)]

# Function to plot a vector and a matrix side by side
function plot_pair(vec, mat)
    p1 = plot(vec, title="Vector Plot")
    p2 = heatmap(mat, title="Matrix Heatmap")
    plot(p1, p2, layout=(1, 2))
end

using FFMPEG

video_path = "tst/output_video.mp4"  # Path to save the video
fps = 24  # Frames per second

# Create a video
@info "Creating video..."
anim = @animate for i in 1:length(vec_of_vecs)-1
    plot_pair(vec_of_vecs[i], vec_of_mats[i])
end

# Compile the frames into a video
mp4(anim, video_path, fps=fps)
@info "Video created: $video_path"
=#

#=
using Plots
using FFMPEG

# Parameters for the Ornstein-Uhlenbeck process
θ = 3.0  # mean reversion rate
σ = 0.05  # volatility

# Grid settings
x_min, x_max = -4.0, 4.0  # range of x (spatial domain)
dx = 0.01  # spatial step size
dt = 0.01  # time step size
t_max = 5.0  # max time
x = x_min:dx:x_max
t = 0:dt:t_max

# Initialize probability density
p = zeros(length(x), length(t))
p[:, 1] = exp.(-x.^2 / (2*σ^2)) / sqrt(2*pi*σ^2)  # initial condition: standard normal distribution

courant = θ*dt/dx


# Solve the Fokker-Planck equation using finite difference method
for j in 1:length(t)-1
    for i in 2:length(x)-1
        p[i, j+1] = p[i, j]
        p[i, j+1]+= dt * θ * p[i, j]
        p[i, j+1]+= dt * θ * x[i] * (p[i+1, j] - p[i-1, j]) / (2*dx)
        p[i, j+1]+= dt * σ^2 / 2 * (p[i+1, j] - 2*p[i, j] + p[i-1, j]) / dx^2
    end
    global courant = θ*dt/dx
end

# Create animation
anim = @animate for j in 1:length(t)
    plot(x, p[:, j], title="p1,t=$(round(t[j], digits=2)),cn=$courant", xlabel="x", ylabel="Probability Density", legend=false)
end

# Save the animation as a video file
gif(anim, "probability_density_evolution.gif", fps = 24)
=#

using PyCall
using FEniCS

# Create mesh and define function space
nx, ny = 10, 10
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1] + t", t=0.0, degree=2)
bc = DirichletBC(V, u_D, "on_boundary")

# Define initial value
u_n = interpolate(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)

dt = 0.1  # time step
u_ = Function(V)
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
t = 0.0
T = 2.0  # final time
while t < T
    t += dt
    u_D.t = t
    solve(a == L, u_, bc)

    # Update previous solution
    u_n.assign(u_)

    # Optionally, save or visualize the solution here
end

# Final solution
println("Final Solution at t = ", T, ":")
println(u_)
