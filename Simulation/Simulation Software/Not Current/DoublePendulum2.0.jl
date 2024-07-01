#################### Double Pendulum ####################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_fn!(du,u,p,t)
    m₁,m₂,Lₐ,L₁,L₂,c₁,c₂,g = p
    τ₁ =  40 * sin(t*pi/1)
    τ₂ =  80 * sin(t*pi/.125)
    A = [m₁*L₁^2 + m₂*(Lₐ^2 + L₂^2 + 2*Lₐ*L₂*cos(u[2]))     m₂*(L₂^2 + Lₐ*L₂*cos(u[2]));
         m₂*(L₂^2 + Lₐ*L₂*cos(u[2]))                        m₂*L₂^2                       ];
    F = [τ₁ - c₁*u[3] - g*m₂*L₂*sin(u[2]+u[1]) + m₂*Lₐ*L₂*u[4]*(u[4]+2*u[3])*sin(u[2]) - g*(m₁*L₁ + m₂*Lₐ)*sin(u[1]);
         τ₂ - c₂*u[4] - g*m₂*L₂*sin(u[2]+u[1]) - m₂*Lₐ*L₂*u[3]^2*sin(u[2])                                              ]
    B = inv(A)*F;

    du[1:2] = u[3:4]
    du[3:4] = B[1:2]

    SA[du]
end

### Convert to Cartesian ###
function X_Y(sol,p)
    n =  length(sol)
    u = zeros(2,n)
    for i in 1:n
        u[:,i] = sol.u[i,1][1:2]
    end

     # julia doc string
     Lₐ,Lᵦ,L₁,L₂ = p
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     n = length(u₁)

     x = zeros(5,n)
     x[2,:] = L₁*sin.(u₁)
     x[3,:] = Lₐ*sin.(u₁)
     x[4,:] = x[3,:] + L₂*sin.(u₂)
     x[5,:] = x[3,:] + Lᵦ*sin.(u₂)

     y = zeros(5,n)
     y[2,:] = -L₁*cos.(u₁)
     y[3,:] = -Lₐ*cos.(u₁)
     y[4,:] = y[3,:] + -L₂*cos.(u₂)
     y[5,:] = y[3,:] + -Lᵦ*cos.(u₂)
     x,y
end


##################### Start of Code #####################
### Initial Conditions ###
x₀ = [1, -2, 0, 0];
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,g = (8,16,4,5,3,3,7,3,9.81) 
p = (m₁,m₂,Lₐ,L₁,L₂,c₁,c₂,g)

####### ODE Solver #######
tspan = (0.0, 100.0);
prob = ODEProblem(ode_fn!,x₀,tspan,p);
sol = solve(prob,Tsit5());

##### Interpret Data #####
p = Lₐ,Lᵦ,L₁,L₂
x,y = X_Y(sol,p)

##### Phase Plot #####
plot(x, y)

####### Create GIF #######
n =  length(sol)
@gif for i ∈ 1:n
     plot(x[:,i],y[:,i],ylim = [-10,0],xlim=[-5,5],xlabel="X Position", ylabel="Y Position", 
        markercolors = [:black, :red, :black, :red, :blue], shape = [:hline,:circle, :cross, :circle, :circle],
        markersize = [300,8,6,8,16], linewidth = 4)
end every 1
