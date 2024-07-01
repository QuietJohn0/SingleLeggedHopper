################## Double Pendulum 3.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_fn!(du,u,p,t)
     τ =  [0 0]'
     M = [p.Lᵦ^2*p.m₁ + 2*p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂ + (p.L₂-p.Lₐ)^2*p.m₁        p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂ + p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]);
          p.Lᵦ^2*p.m₁ +   p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂                             p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂];
    
     V = [p.c₁*u[5] - u[6]*(2*u[5] + u[6])*p.Lᵦ*p.m₁*(p.Lₐ-p.L₂)*sin(u[2]) - p.m₁*p.g*(p.Lᵦ*sin(u[1] + u[2]) + (p.Lₐ-p.L₂)*sin(u[1])) - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2]);
          p.c₂*u[6] + u[5]^2*p.m₁*p.Lᵦ*(p.Lₐ-p.L₂)*sin(u[2])               - p.m₁*p.g*p.Lᵦ*sin(u[1] + u[2])                           - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2])]
                                                                            
     C = [1 0;
          0 1;
          -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2]);
          -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])]
   
     D = [0
          0
          p.Lᵦ*(u[5]+u[6])^2*sin(u[1]+u[2]) + p.Lₐ*u[5]^2*sin(u[1]);
         -p.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - p.Lₐ*u[5]^2*cos(u[1])]

     q̈ = C*(M\(τ-V)) + D

     du[1:4] = u[5:8]
     du[5:8] = q̈
 
     SA[du]
end

### Convert to Cartesian ###
function X_Y(u,p)

     Lₐ,Lᵦ,L₁,L₂ = p
     xₒ = u[3,:]
     yₒ = u[4,:]
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     x = zeros(5,n)
     x[1,:] = xₒ
     x[2,:] = x[1,:] + L₁*sin.(u₁)
     x[3,:] = x[1,:] + Lₐ*sin.(u₁)
     x[4,:] = x[3,:] + L₂*sin.(u₂)
     x[5,:] = x[3,:] + Lᵦ*sin.(u₂)

     y = zeros(5,n)
     y[1,:] = yₒ
     y[2,:] = y[1,:] + -L₁*cos.(u₁)
     y[3,:] = y[1,:] + -Lₐ*cos.(u₁)
     y[4,:] = y[3,:] + -L₂*cos.(u₂)
     y[5,:] = y[3,:] + -Lᵦ*cos.(u₂)
     x,y
end
#=function X_Y(u,p)

     Lₐ,Lᵦ,L₁,L₂ = p
     xₒ = 0
     yₒ = 0
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     x = zeros(5,n)
     x[5,:] = xₒ * ones(1,n)
     x[4,:] = x[5,:] - (Lᵦ-L₂)*sin.(u₂)
     x[3,:] = x[5,:] -      Lᵦ*sin.(u₂)
     x[2,:] = x[3,:] - (Lₐ-L₁)*sin.(u₁)
     x[1,:] = x[3,:] -      L₁*sin.(u₁)

     y = zeros(5,n)
     y[1,:] = yₒ * ones(1,n)
     y[4,:] = y[5,:] + (Lᵦ-L₂)*cos.(u₂)
     y[3,:] = y[5,:] +      Lᵦ*cos.(u₂)
     y[2,:] = y[3,:] + (Lₐ-L₁)*cos.(u₁)
     y[1,:] = y[3,:] +      L₁*cos.(u₁)
     x,y
end=#


##################### Start of Code #####################
### Initial Conditions ###
x₀ = [0, .4, -1.8, 8, 0, 0, 0, 0];
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,cₓ,c🥰,g = (8,16,4,5,3,3,20,20,1,2, 1.81) 
p = (m₁=m₁,m₂=m₂,Lₐ=Lₐ,Lᵦ=Lᵦ,L₁=L₁,L₂=L₂,c₁=c₁,c₂=c₂,cₓ=cₓ,c🥰=c🥰, g=g)

####### ODE Solver #######
tspan = (0.0, 15.0);
prob = ODEProblem(ode_fn!,x₀,tspan,p);
sol = solve(prob,AutoTsit5(Rosenbrock23()),dt=1);

##### Interpret Data #####
n =  length(sol)
dim = length(sol.u[1])
u = zeros(dim,n)
for i in 1:n
    u[:,i] = sol.u[i,1][1:dim]
end

p = Lₐ,Lᵦ,L₁,L₂
x,y = X_Y(u,p)

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

####### Create GIF #######
n =  length(sol)
@gif for i ∈ 1:n
     hi = sol.u[i][4] - Lₐ*cos(sol.u[i][1]) - Lᵦ*cos(sol.u[i][1] + sol.u[i][2])
     plot(x[:,i],y[:,i],ylim = [-10,10],xlim=[-10,10],xlabel="X Position", ylabel="Y Position", 
        markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
        markersize = [3,8,6,8,16], linewidth = 4, title = "Single Leg In Stance", aspect_ratio=:equal)
end every 1
