################## Double Pendulum 3.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_S_fn!(du,u,p,t)
     τ =  [Controller_S(u,p,t); 0]

     v1 = 2*p.Lᵦ*p.L₂*p.m₁*u[4]*u[5]*sin(u[2]) + p.Lᵦ*p.L₂*p.m₁*u[5]^2*sin(u[2]) - 2*p.Lᵦ*p.Lₐ*p.m₁*u[4]*u[5]*sin(u[2]) - p.Lᵦ*p.Lₐ*p.m₁*u[5]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₂*p.g*p.m₁*sin(u[1]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.g*p.m₁*sin(u[1]) + p.c₁*u[4]
     v2 = -p.Lᵦ*p.L₂*p.m₁*u[4]^2*sin(u[2]) + p.Lᵦ*p.Lₐ*p.m₁*u[4]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.c₂*u[5]
     v3 = -p.Lᵦ*p.m₁*u[4]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₁*u[4]*u[5]*cos(u[1] + u[2]) - p.Lᵦ*p.m₁*u[5]^2*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[4]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₂*u[4]*u[5]*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[5]^2*cos(u[1] + u[2]) + p.L₂*p.m₁*u[4]^2*cos(u[1]) + p.L₂*p.m₂*u[4]^2*cos(u[1] + u[2]) + 2*p.L₂*p.m₂*u[4]*u[5]*cos(u[1] + u[2]) + p.L₂*p.m₂*u[5]^2*cos(u[1] + u[2]) - p.Lₐ*p.m₁*u[4]^2*cos(u[1]) + p.g*p.m₁ + p.g*p.m₂ + p.kₛ*u[3]
     
     m11 = p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - 2*p.Lᵦ*p.L₂*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + 2*p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₂^2*p.m₁ + p.L₂^2*p.m₂ - 2*p.L₂*p.Lₐ*p.m₁ + p.Lₐ^2*p.m₁
     m12 = p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - p.Lᵦ*p.L₂*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₂^2*p.m₂
     m13 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.m₁*sin(u[1])
     m21 = p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - p.Lᵦ*p.L₂*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₂^2*p.m₂
     m22 = p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - 2*p.Lᵦ*p.L₂*p.m₂ + p.L₂^2*p.m₂
     m23 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₂*sin(u[1] + u[2])
     m31 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.m₁*sin(u[1])
     m32 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₂*sin(u[1] + u[2])
     m33 = p.m₁ + p.m₂


     M = [m11 m12 m13;
          m21 m22 m23;
          m31 m32 m33]
    
     V = [v1;
          v2;
          v3]

     q̈ = M\(τ-V)

     du[1:3] = u[4:6]
     du[4:6] = q̈
 
     SA[du]
end

### Convert to Cartesian ###
function X_Y(u,p)
     xₒ = 0
     yₒ = u[3,:]
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     x = zeros(5,n)
     x[5,:] = xₒ * ones(1,n)
     x[4,:] = x[5,:] - (p.Lᵦ-p.L₂)*sin.(u₂)
     x[3,:] = x[5,:] -       p.Lᵦ *sin.(u₂)
     x[2,:] = x[3,:] - (p.Lₐ-p.L₁)*sin.(u₁)
     x[1,:] = x[3,:] -       p.L₁ *sin.(u₁)

     y = zeros(5,n)
     y[5,:] = yₒ
     y[4,:] = y[5,:] + (p.Lᵦ-p.L₂)*cos.(u₂)
     y[3,:] = y[5,:] +        p.Lᵦ*cos.(u₂)
     y[2,:] = y[3,:] + (p.Lₐ-p.L₁)*cos.(u₁)
     y[1,:] = y[3,:] +        p.L₁*cos.(u₁)
     x,y
end

function Controller_S(u,p,t)
     # Force
     Fd = p.A*sin(p.w*t + p.ϕ)
     Fa = p.kₛ*u[3]
     e_F = Fd - Fa

     # Velecity at CM
     Vd = p.Vcmd
     Va = Vcm(p,u)
     e_V = Vd - Va

     e = [e_F;e_V]

     # Force derivative
     Ḟd = p.A*p.w*cos(p.w*t + p.ϕ)
     Ḟa = p.kₛ*u[6]
     ė_F = Ḟd - Ḟa

     # Acceteration at CM
     V̇d = 0
     V̇a = Acm(p,u)
     ė_V = V̇d - V̇a

     ė = [ė_F;ė_V]

     τ = p.Kp_s*e + p.Kd_s*ė
end

function Vcm(p,u)
     u = Aug_S(u,0)
     J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
          ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
     Pcm = J*u[5:8]
     ẋcm = Pcm[1]
end
function Acm(p,u)

     u = Aug_S(u,0)
     J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
          ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
     
     Vⱼ= [-u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])) + -u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*sin(u[1]+u[2]);
           u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])) +  u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*cos(u[1]+u[2])]/(p.m₁ + p.m₂)
     
     Pcm = J*u[5:8] + Vⱼ
     ẋcm = Pcm[1]
end

function Aug_S(u, x₀)
     C = [1 0 0;
          0 1 0;
         -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2])    0;
         -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])    1;]
       
     posO = [x₀ + -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])
                   p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1]) + u[3]]
 
     u = [u[1:2]; posO; C*u[4:6]; u[3]; u[6]]
end

mutable struct System_C
     m₁::Float64
     m₂::Float64
     Lₐ::Float64
     Lᵦ::Float64
     L₁::Float64
     L₂::Float64
     c₁::Float64
     c₂::Float64
     kₛ::Float64
     g::Float64

     C::Matrix{Float64}
     tf::Float64
     Kp_f::Matrix{Float64}
     Kd_f::Matrix{Float64}

     A::Float64
     w::Float64
     ϕ::Float64
     Vcmd::Float64
     Kp_s::Matrix{Float64}
     Kd_s::Matrix{Float64}
end

##################### Start of Code #####################
### Initial Conditions ###
x₀ = [1, -2, 0, 0, 0, 0];

### Create System ###
p = System_C(16,8,.2,.18,.05,.06,20,20,5E6,9.81,zeros(4, 2),1,zeros(2, 2),zeros(2, 2),50,50*2*pi,0,1.5,zeros(2, 2),zeros(2, 2))

### Set up Controller ###
p.A = 50
p.w = 50*2*pi
p.ϕ = 0
p.Vcmd = 1.5
p.Kp_s = [.1 .2; .2 .3]
p.Kd_s = [.5 .6; .1 .1]

####### ODE Solver #######
tspan = (0.0, 5);
prob = ODEProblem(ode_S_fn!,x₀,tspan,p)
condition(u,t,integrator) = u[3]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb, abstol=1e-10,reltol=1e-10, saveat = .0001)

##### Interpret Data #####
n =  length(sol)
dim = length(sol.u[1])
u = zeros(dim,n)
for i in 1:n
    u[:,i] = sol.u[i,1][1:dim]
end

x,y = X_Y(u,p)

open("test2.txt", "w") do file
     write(file, string(y))
end

############### Phase Plot ###############
##### X V.S Y #####
plt = plot(x[1,:], y[1,:], lc=:black, lw=2)
plot!(x[3,:], y[3,:], lc=:red, lw=2)
plot!(x[5,:], y[5,:], lc=:blue, lw=2)
plot!([x[1,1]],[y[1,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[3,1]],[y[3,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[5,1]],[y[5,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Y and X Position of Orgin")
xlabel!("x")
ylabel!("y")
display(plt)

##### θ₁ V.S θ̇₁ #####
plt = plot(u[1,:], u[4,:], lc=:black, lw=2)
plot!([u[1,1]],[u[4,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₁ V.S θ̇₁")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)

##### θ₃ V.S θ̇₃ #####
plt = plot(u[2,:], u[5,:], lc=:black, lw=2)
plot!([u[2,1]],[u[5,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₃ V.S θ̇₃")
xlabel!("θ₁")
ylabel!("θ̇₃")
display(plt)

##### yₛ V.S ẏₛ #####
plt = plot(u[3,:], u[6,:], lc=:black, lw=2)
plot!([u[3,1]],[u[6,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot yₛ V.S ẏₛ")
xlabel!("yₛ")
ylabel!("ẏₛ")
display(plt)

##### θ₁ V.S θ̇₁ #####


# Force
Fd = p.A*sin.(p.w*sol.t)

plt = plot(sol.t, u[3,:]*p.kₛ, lc=:black, lw=2)
plt = plot!(sol.t, Fd, lc=:blue, lw=2)
plot!([sol.t[1]],[u[3,1]*p.kₛ],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot t V.S F")
xlabel!("t")
ylabel!("F")
display(plt)

####### Create GIF #######
n =  length(sol)
s = p
@gif for i ∈ 1:n
     #hi = sol.u[i][4] - p.Lₐ*cos(sol.u[i][1]) - p.Lᵦ*cos(sol.u[i][1]+ sol.u[i][2])
     hi = sol.u[i][4] - p.Lₐ*cos(sol.u[i][1]) - p.Lᵦ*cos(sol.u[i][1]+sol.u[i][2])
     plot(x[:,i],y[:,i],ylim = [-.5,.5],xlim=[-.8,.5],xlabel="X Position (m)", ylabel="Y Position (m)", 
        markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
        markersize = [3,8,6,8,16], linewidth = 4, title = sol.u[i][3], aspect_ratio=:equal)
 end every 10
