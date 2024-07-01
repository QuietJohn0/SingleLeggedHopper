################## Double Pendulum 3.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_F_fn!(du,u,p,t)
    τ =  [Controller_F(u,p,t); 0; 0]

    m11 = p.L₁^2*p.m₁ + p.L₂^2*p.m₂ + 2*p.L₂*p.Lₐ*p.m₂*cos(u[2]) + p.Lₐ^2*p.m₂
    m12 = p.L₂^2*p.m₂ + p.L₂*p.Lₐ*p.m₂*cos(u[2])
    m13 = p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) + p.Lₐ*p.m₂*cos(u[1])
    m14 = p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.m₂*sin(u[1])
    m21 = p.L₂^2*p.m₂ + p.L₂*p.Lₐ*p.m₂*cos(u[2])
    m22 = p.L₂^2*p.m₂
    m23 = p.L₂*p.m₂*cos(u[1] + u[2])
    m24 = p.L₂*p.m₂*sin(u[1] + u[2])
    m31 = p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) + p.Lₐ*p.m₂*cos(u[1])
    m32 = p.L₂*p.m₂*cos(u[1] + u[2])
    m33 = p.m₁ + p.m₂
    m34 = 0
    m41 = p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.m₂*sin(u[1])
    m42 = p.L₂*p.m₂*sin(u[1] + u[2])
    m43 = 0
    m44 = p.m₁ + p.m₂

    v1 = p.L₁*p.g*p.m₁*sin(u[1]) - 2*p.L₂*p.Lₐ*p.m₂*u[5]*u[6]*sin(u[2]) - p.L₂*p.Lₐ*p.m₂*u[6]^2*sin(u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.g*p.m₂*sin(u[1]) + p.c₁*u[5]
    v2 = p.L₂*p.Lₐ*p.m₂*u[5]^2*sin(u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.c₂*u[6]
    v3 = -p.L₁*p.m₁*u[5]^2*sin(u[1]) - p.L₂*p.m₂*u[5]^2*sin(u[1] + u[2]) - 2*p.L₂*p.m₂*u[5]*u[6]*sin(u[1] + u[2]) - p.L₂*p.m₂*u[6]^2*sin(u[1] + u[2]) - p.Lₐ*p.m₂*u[5]^2*sin(u[1])
    v4 = p.L₁*p.m₁*u[5]^2*cos(u[1]) + p.L₂*p.m₂*u[5]^2*cos(u[1] + u[2]) + 2*p.L₂*p.m₂*u[5]*u[6]*cos(u[1] + u[2]) + p.L₂*p.m₂*u[6]^2*cos(u[1] + u[2]) + p.Lₐ*p.m₂*u[5]^2*cos(u[1]) + p.g*p.m₁ + p.g*p.m₂



    M = [m11 m12 m13 m14;
         m21 m22 m23 m24;
         m31 m32 m33 m34;
         m41 m42 m43 m44]
    
    V = [v1;
         v2;
         v3;
         v4]

    q̈ = M\(τ-V)

    du[1:4] = u[5:8]
    du[5:8] = q̈

    SA[du]
end

### Convert to Cartesian ###
function X_Y(u,p)
     xₒ = u[3,:]
     yₒ = u[4,:]
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     x = zeros(5,n)
     x[1,:] = xₒ
     x[2,:] = x[1,:] + p.L₁*sin.(u₁)
     x[3,:] = x[1,:] + p.Lₐ*sin.(u₁)
     x[4,:] = x[3,:] + p.L₂*sin.(u₂)
     x[5,:] = x[3,:] + p.Lᵦ*sin.(u₂)

     y = zeros(5,n)
     y[1,:] = yₒ
     y[2,:] = y[1,:] + -p.L₁*cos.(u₁)
     y[3,:] = y[1,:] + -p.Lₐ*cos.(u₁)
     y[4,:] = y[3,:] + -p.L₂*cos.(u₂)
     y[5,:] = y[3,:] + -p.Lᵦ*cos.(u₂)
     x,y
end
### Controller ###
function Controller_F(u,p,t)
     θd = Traj(p,t)

     
     # Position θ₁
     θ₁d = θd[1,1]
     θ₁a = u[1]
     e_θ₁ = θ₁d - θ₁a

     # Position θ₃
     θ₃d = θd[1,2]
     θ₃a = u[2]
     e_θ₃ = θ₃d - θ₃a

     e = [e_θ₁;e_θ₃]

     # Position Derivative, θ̇₁
     θ̇₁d = θd[2,1]
     θ̇₁a = u[5]
     ė_θ₁ = θ̇₁d - θ̇₁a

     # Acceteration at CM, θ̇₃
     θ̇₃d = θd[2,2]
     θ̇₃a = u[6]
     ė_θ₃ = θ̇₃d - θ̇₃a

     ė = [ė_θ₁;ė_θ₃]

     τ = p.Kp_f*e + p.Kd_f*ė
end
function Traj(p,t)
     if p.tf < t
          t = p.tf
     end

     C = p.C
     θd = [1 t t^2 t^3;
           0 1 2*t 3*t^2]*C
end
function TrajPlan(u,tf, setpoint)
     M_t =[1   0    0    0;
           0   1    0    0;
           1   tf   tf^2 tf^3;
           0   1    2*tf 3*tf^2]
     M_IC = [u[1:2]'; u[5:6]'; setpoint]
     C = M_t\M_IC
end

### System Structure ###
mutable struct System
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
x₀ = [0, 0, 0, 5, .8, -.9, 0, 0]
x₀ = [-.9, .9, 0, 65, -.3, .8, -.5, 0]

### Create System ###
p = System(8,16,4,5,3,3,20,20,5E6,9.81,zeros(4, 2),1,zeros(2, 2),zeros(2, 2),50,50*2*pi,0,1.5,zeros(2, 2),zeros(2, 2))

### Set up Fight Controllers ###
setpoint = [0 0; 0 0]
p.tf = 3
p.C = TrajPlan(x₀,p.tf, setpoint)
p.Kp_f = [650 400; 300 300]
p.Kd_f = [50 60; 10 10]

####### ODE Solver #######
tspan = (0,5);
prob = ODEProblem(ode_F_fn!,x₀,tspan,p);
condition(u,t,integrator) = u[4] - integrator.p.Lₐ*cos(u[1]) - integrator.p.Lᵦ*cos(u[1]+u[2])
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb, abstol=1e-10,reltol=1e-10, saveat = .01);

##### Interpret Data #####
n =  length(sol)
dim = length(sol.u[1])
u = zeros(dim,n)
for i in 1:n
    u[:,i] = sol.u[i,1][1:dim]
end

x,y = X_Y(u,p)
t = sol.t

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
idx = Int(dim/2+1)
plt = plot(u[1,:], u[idx,:], lc=:black, lw=2)
plot!([u[1,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₁ V.S θ̇₁")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)

##### θ₃ V.S θ̇₃ #####
idx = Int(dim/2+2)
plt = plot(u[2,:], u[idx,:], lc=:black, lw=2)
plot!([u[2,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₃ V.S θ̇₃")
xlabel!("θ₁")
ylabel!("θ̇₃")
display(plt)

##### x V.S ẋ #####
idx = Int(dim/2+3)
plt = plot(u[3,:], u[idx,:], lc=:black, lw=2)
plot!([u[3,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot x V.S ẋ")
xlabel!("x")
ylabel!("ẋ")
display(plt)

##### y V.S ẏ #####
idx = Int(dim/2+4)
plt = plot(u[4,:], u[idx,:], lc=:black, lw=2)
plot!([u[4,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot y V.S ẏ")
xlabel!("y")
ylabel!("ẏ")
display(plt)

##### θ₁,θ₃  V.S t #####

θ = [ones(n,1) t t.^2 t.^3;]*p.C

plt = plot(t,u[1,:], lc=:black, lw=2)
plot!(t, θ[:,1],markercolors = :red, markersize = 3)
plot!(t, u[2,:], lc=:black, lw=2)
plot!(t, θ[:,2],markercolors = :red, markersize = 3)
title!("Phase Plot t V.S θ")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)


####### Create GIF #######
n =  length(sol)
title = "Single Leg In Flight"
@gif for i ∈ 1:n
     plot(x[:,i],y[:,i],ylim = [-2,80],xlim=[-15,15],xlabel="X Position", ylabel="Y Position", 
        markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
        markersize = [3,8,6,8,16], linewidth = 4, title = sol.t[i],aspect_ratio=:equal)
end every 5
