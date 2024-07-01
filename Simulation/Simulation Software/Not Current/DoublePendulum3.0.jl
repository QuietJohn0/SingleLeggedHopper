################## Double Pendulum 3.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_fn!(du,u,p,t)
     τ =  [0 0 0 0]'
     
     M = [p.m₁*p.L₁^2 + p.m₂*(p.Lₐ^2 + p.L₂^2 + 2*p.Lₐ*p.L₂*cos(u[2]))      p.m₂*(p.L₂^2 + p.Lₐ*p.L₂*cos(u[2]))      (p.m₁*p.L₁+p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])    (p.m₁*p.L₁+p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]);
         p.m₂*(p.L₂^2 + p.Lₐ*p.L₂*cos(u[2]))                               p.m₂*p.L₂^2                               p.m₂*p.L₂*cos(u[1]+u[2])                                    p.m₂*p.L₂*sin(u[1]+u[2]);
         (p.m₁*p.L₁+p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])          p.m₂*p.L₂*cos(u[1]+u[2])                   p.m₁ + p.m₂                                                0;
         (p.m₁*p.L₁+p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])          p.m₂*p.L₂*sin(u[1]+u[2])                   0                                                          p.m₁ + p.m₂];
     V = [p.c₁*u[5] - p.m₂*p.Lₐ*p.L₂*u[6]*(u[6]+2*u[5])*sin(u[2]) + p.g*p.m₂*p.L₂*sin(u[1]+u[2]) + p.g*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]);
          p.c₂*u[6] + p.m₂*p.Lₐ*p.L₂*u[5]^2*sin(u[2])             + p.g*p.m₂*p.L₂*sin(u[1]+u[2]);
         -u[5]^2*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[2])      - (u[5]+u[6])^2*p.m₂*p.L₂*sin(u[1]+u[2]);
          u[5]^2*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[2])      + (u[5]+u[6])^2*p.m₂*p.L₂*cos(u[1]+u[2]) + p.g*(p.m₁ + p.m₂)];
     
     q̈ = M\(τ-V)

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


##################### Start of Code #####################
### Initial Conditions ###
x₀ = [-.9, .9, 0, 0, -.3, .8, -.5, 7];
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,g = (8,16,4,5,3,3,20,20, 1.81) 
p = (m₁=m₁,m₂=m₂,Lₐ=Lₐ,L₁=L₁,L₂=L₂,c₁=c₁,c₂=c₂, g=g)

####### ODE Solver #######
tspan = (0.0, 10.0);
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


####### Create GIF #######
n =  length(sol)
@gif for i ∈ 1:n
     plot(x[:,i],y[:,i],ylim = [-10,20],xlim=[-15,15],xlabel="X Position", ylabel="Y Position", 
        markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
        markersize = [2,2,1,2,4], linewidth = 1, title = "Single Leg In Flight",aspect_ratio=:equal)
end every 1
