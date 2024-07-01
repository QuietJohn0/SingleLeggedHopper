################## Double Pendulum Combined 1.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_F_fn!(du,u,p,t)
    τ =  [0 0 0 0]'
    s = p

    m11 = s.L₁^2*s.m₁ + s.L₂^2*s.m₂ + 2*s.L₂*s.Lₐ*s.m₂*cos(u[2]) + s.Lₐ^2*s.m₂
    m12 = s.L₂^2*s.m₂ + s.L₂*s.Lₐ*s.m₂*cos(u[2])
    m13 = s.L₁*s.m₁*cos(u[1]) + s.L₂*s.m₂*cos(u[1] + u[2]) + s.Lₐ*s.m₂*cos(u[1])
    m14 = s.L₁*s.m₁*sin(u[1]) + s.L₂*s.m₂*sin(u[1] + u[2]) + s.Lₐ*s.m₂*sin(u[1])
    m21 = s.L₂^2*s.m₂ + s.L₂*s.Lₐ*s.m₂*cos(u[2])
    m22 = s.L₂^2*s.m₂
    m23 = s.L₂*s.m₂*cos(u[1] + u[2])
    m24 = s.L₂*s.m₂*sin(u[1] + u[2])
    m31 = s.L₁*s.m₁*cos(u[1]) + s.L₂*s.m₂*cos(u[1] + u[2]) + s.Lₐ*s.m₂*cos(u[1])
    m32 = s.L₂*s.m₂*cos(u[1] + u[2])
    m33 = s.m₁ + s.m₂
    m34 = 0
    m41 = s.L₁*s.m₁*sin(u[1]) + s.L₂*s.m₂*sin(u[1] + u[2]) + s.Lₐ*s.m₂*sin(u[1])
    m42 = s.L₂*s.m₂*sin(u[1] + u[2])
    m43 = 0
    m44 = s.m₁ + s.m₂

    v1 = s.L₁*s.g*s.m₁*sin(u[1]) - 2*s.L₂*s.Lₐ*s.m₂*u[5]*u[6]*sin(u[2]) - s.L₂*s.Lₐ*s.m₂*u[6]^2*sin(u[2]) + s.L₂*s.g*s.m₂*sin(u[1] + u[2]) + s.Lₐ*s.g*s.m₂*sin(u[1]) + s.c₁*u[5]
    v2 = s.L₂*s.Lₐ*s.m₂*u[5]^2*sin(u[2]) + s.L₂*s.g*s.m₂*sin(u[1] + u[2]) + s.c₂*u[6]
    v3 = -s.L₁*s.m₁*u[5]^2*sin(u[1]) - s.L₂*s.m₂*u[5]^2*sin(u[1] + u[2]) - 2*s.L₂*s.m₂*u[5]*u[6]*sin(u[1] + u[2]) - s.L₂*s.m₂*u[6]^2*sin(u[1] + u[2]) - s.Lₐ*s.m₂*u[5]^2*sin(u[1])
    v4 = s.L₁*s.m₁*u[5]^2*cos(u[1]) + s.L₂*s.m₂*u[5]^2*cos(u[1] + u[2]) + 2*s.L₂*s.m₂*u[5]*u[6]*cos(u[1] + u[2]) + s.L₂*s.m₂*u[6]^2*cos(u[1] + u[2]) + s.Lₐ*s.m₂*u[5]^2*cos(u[1]) + s.g*s.m₁ + s.g*s.m₂



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
function ode_S_fn!(du,u,p,t)
    τ =  [80 -150 0]'
    s = p

    m11 = s.Lᵦ^2*s.m₁ + s.Lᵦ^2*s.m₂ - 2*s.Lᵦ*s.L₂*s.m₁*cos(u[2]) - 2*s.Lᵦ*s.L₂*s.m₂ + 2*s.Lᵦ*s.Lₐ*s.m₁*cos(u[2]) + s.L₂^2*s.m₁ + s.L₂^2*s.m₂ - 2*s.L₂*s.Lₐ*s.m₁ + s.Lₐ^2*s.m₁
    m12 = s.Lᵦ^2*s.m₁ + s.Lᵦ^2*s.m₂ - s.Lᵦ*s.L₂*s.m₁*cos(u[2]) - 2*s.Lᵦ*s.L₂*s.m₂ + s.Lᵦ*s.Lₐ*s.m₁*cos(u[2]) + s.L₂^2*s.m₂
    m13 = -s.Lᵦ*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.m₂*sin(u[1] + u[2]) + s.L₂*s.m₁*sin(u[1]) + s.L₂*s.m₂*sin(u[1] + u[2]) - s.Lₐ*s.m₁*sin(u[1])
    m21 = s.Lᵦ^2*s.m₁ + s.Lᵦ^2*s.m₂ - s.Lᵦ*s.L₂*s.m₁*cos(u[2]) - 2*s.Lᵦ*s.L₂*s.m₂ + s.Lᵦ*s.Lₐ*s.m₁*cos(u[2]) + s.L₂^2*s.m₂
    m22 = s.Lᵦ^2*s.m₁ + s.Lᵦ^2*s.m₂ - 2*s.Lᵦ*s.L₂*s.m₂ + s.L₂^2*s.m₂
    m23 = -s.Lᵦ*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.m₂*sin(u[1] + u[2]) + s.L₂*s.m₂*sin(u[1] + u[2])
    m31 = -s.Lᵦ*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.m₂*sin(u[1] + u[2]) + s.L₂*s.m₁*sin(u[1]) + s.L₂*s.m₂*sin(u[1] + u[2]) - s.Lₐ*s.m₁*sin(u[1])
    m32 = -s.Lᵦ*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.m₂*sin(u[1] + u[2]) + s.L₂*s.m₂*sin(u[1] + u[2])
    m33 = s.m₁ + s.m₂
    
    v1 = 2*s.Lᵦ*s.L₂*s.m₁*u[4]*u[5]*sin(u[2]) + s.Lᵦ*s.L₂*s.m₁*u[5]^2*sin(u[2]) - 2*s.Lᵦ*s.Lₐ*s.m₁*u[4]*u[5]*sin(u[2]) - s.Lᵦ*s.Lₐ*s.m₁*u[5]^2*sin(u[2]) - s.Lᵦ*s.g*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.g*s.m₂*sin(u[1] + u[2]) + s.L₂*s.g*s.m₁*sin(u[1]) + s.L₂*s.g*s.m₂*sin(u[1] + u[2]) - s.Lₐ*s.g*s.m₁*sin(u[1]) + s.c₁*u[4]
    v2 = -s.Lᵦ*s.L₂*s.m₁*u[4]^2*sin(u[2]) + s.Lᵦ*s.Lₐ*s.m₁*u[4]^2*sin(u[2]) - s.Lᵦ*s.g*s.m₁*sin(u[1] + u[2]) - s.Lᵦ*s.g*s.m₂*sin(u[1] + u[2]) + s.L₂*s.g*s.m₂*sin(u[1] + u[2]) + s.c₂*u[5]
    v3 = -s.Lᵦ*s.m₁*u[4]^2*cos(u[1] + u[2]) - 2*s.Lᵦ*s.m₁*u[4]*u[5]*cos(u[1] + u[2]) - s.Lᵦ*s.m₁*u[5]^2*cos(u[1] + u[2]) - s.Lᵦ*s.m₂*u[4]^2*cos(u[1] + u[2]) - 2*s.Lᵦ*s.m₂*u[4]*u[5]*cos(u[1] + u[2]) - s.Lᵦ*s.m₂*u[5]^2*cos(u[1] + u[2]) + s.L₂*s.m₁*u[4]^2*cos(u[1]) + s.L₂*s.m₂*u[4]^2*cos(u[1] + u[2]) + 2*s.L₂*s.m₂*u[4]*u[5]*cos(u[1] + u[2]) + s.L₂*s.m₂*u[5]^2*cos(u[1] + u[2]) - s.Lₐ*s.m₁*u[4]^2*cos(u[1]) + s.g*s.m₁ + s.g*s.m₂ + s.kₛ*u[3]



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

######## Controllers ########
function PDPosition(u,d)
    Kp = [2 4; 5 7]
    Kd = [2 4; 5 7]
    e = [u[1:4]-d[1:4]]
    τ = [Kp Kd]*e
end
function PDForceVcm(u,d,dt)
    Kp = [2 4; 5 7]
    Kd = [2 4; 5 7]
    e = [u[1:4]-d[1:4]]
    τ = [Kp Kd]*e
end

##### Convert Functions #####
function GRF2Torque(s,u,Bf)
    ## This function Will be used to convert desired GRF to expected torque
    M = [s.Lᵦ^2*s.m₁ + 2*s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂ + (s.L₂-s.Lₐ)^2*s.m₁        s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂ + s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]);
         s.Lᵦ^2*s.m₁ +   s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂                             s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂];
    V = [s.c₁*u[5] - u[6]*(2*u[5] + u[6])*s.Lᵦ*s.m₁*(s.Lₐ-s.L₂)*sin(u[2]) - s.m₁*s.g*(s.Lᵦ*sin(u[1] + u[2]) + (s.Lₐ-s.L₂)*sin(u[1])) - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2]);
         s.c₂*u[6] + u[5]^2*s.m₁*s.Lᵦ*(s.Lₐ-s.L₂)*sin(u[2])               - s.m₁*s.g*s.Lᵦ*sin(u[1] + u[2])                           - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2])]
    A = [-s.Lᵦ*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*cos(u[1])   -s.Lᵦ*cos(u[1]+u[2]);
         -s.Lᵦ*sin(u[1]+u[2]) - (s.Lₐ-s.L₁)*sin(u[1])   -s.Lᵦ*sin(u[1]+u[2]);
         -(s.Lᵦ-s.L₂)*cos(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*cos(u[1]+u[2]) ;
         -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])]
    F = [ s.Lᵦ*(u[5]+u[6])^2*sin(u[1]+u[2]) + (s.Lₐ-s.L₁)*u[5]^2*sin(u[1]);
         -s.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*u[5]^2*cos(u[1]);
          (s.Lᵦ-s.L₂)*(u[5]+u[6])^2*sin(u[1]+u[2])
         -(s.Lᵦ-s.L₂)*(u[5]+u[6])^2*cos(u[1]+u[2])]
    mₘ = [s.m₁ 0    s.m₂ 0;
          0    s.m₁ 0    s.m₂]
    τ = M*inv(mₘ*A)*(Bf - mₘ*F - [0; g*(s.m₁ + s.m₂)]) + V    
end
function Vcm(s,u)
    J = [((s.m₁*s.L₁ + s.m₂*s.Lₐ)*cos(u[1]) + s.m₂*s.L₂*cos(u[1]+u[2]))/(s.m₁ + s.m₂)     s.m₂*s.L₂*cos(u[1]+u[2])/(s.m₁ + s.m₂)    1   0;
         ((s.m₁*s.L₁ + s.m₂*s.Lₐ)*sin(u[1]) + s.m₂*s.L₂*sin(u[1]+u[2]))/(s.m₁ + s.m₂)     s.m₂*s.L₂*sin(u[1]+u[2])/(s.m₁ + s.m₂)    0   1]
    Pcm = J*u[5:8]
    ẋcm = Pcm[1]
end
function impact(u,p)
    s = p
    A = [-s.Lᵦ*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*cos(u[1])   -s.Lᵦ*cos(u[1]+u[2]);
         -s.Lᵦ*sin(u[1]+u[2]) - (s.Lₐ-s.L₁)*sin(u[1])   -s.Lᵦ*sin(u[1]+u[2]);
         -(s.Lᵦ-s.L₂)*cos(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*cos(u[1]+u[2]) ;
         -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])]

    B = [s.L₁*cos(u[1])                       0                     1   0;
         s.L₁*sin(u[1])                       0                     0   1;
         s.Lₐ*cos(u[1]) + s.L₂*cos(u[1]+u[2]) s.L₂*cos(u[1]+u[2])   1   0;
         s.Lₐ*sin(u[1]) + s.L₂*sin(u[1]+u[2]) s.L₂*sin(u[1]+u[2])   0   1]

    M = [-m₁*(s.Lᵦ*cos(u[1]+u[2]) + (s.Lₐ-s.L₁)*cos(u[1]))  -m₁*(s.Lᵦ*sin(u[1]+u[2]) + (s.Lₐ-s.L₁)*sin(u[1]))   -m₂*((s.Lᵦ-s.L₂)*cos(u[1]+u[2]))    -m₂*((s.Lᵦ-s.L₂)*sin(u[1]+u[2]));
         -m₁*(s.Lₐ-s.L₁)*cos(u[1])                          -m₁*(s.Lₐ-s.L₁)*sin(u[1])                            0                                   0]

    C = [1 0;
         0 1;
    -s.Lᵦ*cos(u[1]+u[2]) - s.Lₐ*cos(u[1])  -s.Lᵦ*cos(u[1]+u[2]);
    -s.Lᵦ*sin(u[1]+u[2]) - s.Lₐ*sin(u[1])  -s.Lᵦ*sin(u[1]+u[2])]
    
    u[5:8] = C*inv(M*A)*M*B*u[5:8]
    u
end

##### State Mapping #####
function q_F(u_all, u)
    u_temp = [reshape(u[1],(1,8)) 0 0]
    for i in 2:length(u)
        u_temp = [u_temp; reshape(u[i],(1,8)) 0 0]
    end
    u_all = [u_all; u_temp]
end

function q_S(u_all, u, x₀)
    u_temp = reshape(Aug_S(u[1],x₀),(1,10))
    for i in 2:length(u)
        u_temp = [u_temp; reshape(Aug_S(u[i],x₀),(1,10))]
    end
    u_all = [u_all; u_temp]
end

function Aug_S(u, x₀)
    C = [1 0 0;
         0 1 0;
        -s.Lᵦ*cos(u[1]+u[2]) - s.Lₐ*cos(u[1])  -s.Lᵦ*cos(u[1]+u[2])    0;
        -s.Lᵦ*sin(u[1]+u[2]) - s.Lₐ*sin(u[1])  -s.Lᵦ*sin(u[1]+u[2])    1;]
      
    posO = [x₀ + -s.Lᵦ*sin(u[1]+u[2]) - s.Lₐ*sin(u[1])
                  s.Lᵦ*cos(u[1]+u[2]) + s.Lₐ*cos(u[1]) + u[3]]

    u = [u[1:2]; posO; C*u[4:6]; u[3]; u[6]]
end

### Convert to Cartesian ###
function X_Y(u,s)
    n = length(u[:,1])

    Lₐ,Lᵦ,L₁,L₂ = s
    xₒ = u[:,3]
    yₒ = u[:,4]
    u₁ = u[:,1]
    u₂ = u[:,2] + u[:,1]
    yₛ = u[:,9]

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
# System Parameters #
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,kₛ,g = (16,8,.2,.18,.05,.06,15,15,5E5,9.81)

p = (m₁=m₁,m₂=m₂,Lₐ=Lₐ,Lᵦ=Lᵦ,L₁=L₁,L₂=L₂,c₁=c₁,c₂=c₂,kₛ=kₛ,g=g)

# Controller Stored Values #
tk,d,Vcmd,ωd,Ad,ϕd = (0, [-.4; -.2; 0; 0], .5, .3,50,.2)
c = (tk=tk,d=d,Vcmd=Vcmd,ωd=ωd,Ad=Ad,ϕd=ϕd)

### Initial Conditions ###
x₀ = [-.8 1.6 .24 .5 0 0 0 2];

####### ODE Solver #######
tEnd = 7
affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3]
cb_lift = ContinuousCallback(c_lift,affect!)
c_impact(u,t,integrator) = u[4] - integrator.p.Lₐ*cos(u[1]) - integrator.p.Lᵦ*cos(u[1]+u[2])
cb_impact = ContinuousCallback(c_impact,affect!)

t = 0
t_all = [0]
u_all = [x₀ 0 0]
sw = [0]

#while t <= tEnd
    tspan = (0, tEnd - t)
    prob = ODEProblem(ode_F_fn!,x₀,tspan,p)
    sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = .01)
    
    u_all = q_F(u_all, sol.u)
    t_all = [t_all; sol.t]
    sw = [sw; length(u_all)]

    x₀ = impact(sol.u[end],p)
    t = sol.t[end]

    xᵦ = sol.u[end][3] + p.Lₐ*sin(sol.u[end][1]) + p.Lᵦ*sin(sol.u[end][1]+ sol.u[end][2])

    x₀ = [x₀[1:2]; 0; x₀[5:6]; 0]
    tspan = (0, tEnd - t)
    prob = ODEProblem(ode_S_fn!,x₀,tspan,p)
    sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-10,reltol=1e-10, saveat = .01)

    u_all = q_S(u_all, sol.u, xᵦ)
    t_all = [t_all; sol.t]
    sw = [sw; length(u_all[:,1])]


    x₀ = Aug_S(sol.u[end],xᵦ)[1:8]
#end

display("done")

##### Interpret Data #####
w = Lₐ,Lᵦ,L₁,L₂
x,y = X_Y(u_all,w)

open("test1.txt", "w") do file
    write(file, string(y))
end

############### Phase Plot ###############
##### X V.S Y #####
#=plt = plot(x[1,:], y[1,:], lc=:black, lw=2)
plot!(x[3,:], y[3,:], lc=:red, lw=2)
plot!(x[5,:], y[5,:], lc=:blue, lw=2)
plot!([x[1,1]],[y[1,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[3,1]],[y[3,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[5,1]],[y[5,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Y and X Position of Orgin")
xlabel!("x")
ylabel!("y")
display(plt)=#

##### θ₁ V.S θ̇₁ #####
plt = plot(u_all[:,1], u_all[:,5], lc=:black, lw=2)
plot!([u_all[1,1]],[u_all[1,5]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₁ V.S θ̇₁")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)

##### θ₃ V.S θ̇₃ #####
plt = plot(u_all[:,2], u_all[:,6], lc=:black, lw=2)
plot!([u_all[1,2]],[u_all[1,6]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₃ V.S θ̇₃")
xlabel!("θ₁")
ylabel!("θ̇₃")
display(plt)

##### x V.S ẋ #####
plt = plot(u_all[:,3], u_all[:,7], lc=:black, lw=2)
plot!([u_all[1,3]],[u_all[1,7]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot x V.S ẋ")
xlabel!("x")
ylabel!("ẋ")
display(plt)

##### y V.S ẏ #####
plt = plot(u_all[:,4], u_all[:,8], lc=:black, lw=2)
plot!([u_all[1,4]],[u_all[1,8]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot y V.S ẏ")
xlabel!("y")
ylabel!("ẏ")
display(plt)


####### Create GIF #######
n =  length(u_all[:,1])
s = p
title = "Single Leg With State Transitions"
@gif for i ∈ 1:n
    #hi = sol.u[i][4] - s.Lₐ*cos(sol.u[i][1]) - s.Lᵦ*cos(sol.u[i][1]+ sol.u[i][2])
    plot(x[:,i],y[:,i],ylim = [-.5,.5],xlim=[-.8,.5],xlabel="X Position (m)", ylabel="Y Position (m)", 
       markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
       markersize = [3,8,6,8,16], linewidth = 4, title = u_all[i,9], aspect_ratio=:equal)
    xᵦ = u_all[i,3] + p.Lₐ*sin(u_all[i,1]) + p.Lᵦ*sin(u_all[i,1]+ u_all[i,2])
    plot!([xᵦ], [u_all[i,9]], markercolors = :red, shape = :circle, markersize = 3)

end every 1
