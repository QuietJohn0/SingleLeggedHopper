################## Double Pendulum Combined 4.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays
using XLSX
using DataFrames

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
function ode_S_fn!(du,u,p,t)
    τ =  [Controller_S(u,p,t); 0]

    v1 = 2*p.Lᵦ*p.L₂*p.m₁*u[4]*u[5]*sin(u[2]) + p.Lᵦ*p.L₂*p.m₁*u[5]^2*sin(u[2]) - 2*p.Lᵦ*p.Lₐ*p.m₁*u[4]*u[5]*sin(u[2]) - p.Lᵦ*p.Lₐ*p.m₁*u[5]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₂*p.g*p.m₁*sin(u[1]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.g*p.m₁*sin(u[1]) + p.c₁*u[4]
    v2 = -p.Lᵦ*p.L₂*p.m₁*u[4]^2*sin(u[2]) + p.Lᵦ*p.Lₐ*p.m₁*u[4]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.c₂*u[5]
    v3 = -p.Lᵦ*p.m₁*u[4]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₁*u[4]*u[5]*cos(u[1] + u[2]) - p.Lᵦ*p.m₁*u[5]^2*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[4]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₂*u[4]*u[5]*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[5]^2*cos(u[1] + u[2]) + p.L₂*p.m₁*u[4]^2*cos(u[1]) + p.L₂*p.m₂*u[4]^2*cos(u[1] + u[2]) + 2*p.L₂*p.m₂*u[4]*u[5]*cos(u[1] + u[2]) + p.L₂*p.m₂*u[5]^2*cos(u[1] + u[2]) - p.Lₐ*p.m₁*u[4]^2*cos(u[1]) + p.cₛ*u[6] + p.g*p.m₁ + p.g*p.m₂ + p.kₛ*u[3]
    
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

######## Controllers ########
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

    #  Position Derivative, θ̇₃
    θ̇₃d = θd[2,2]
    θ̇₃a = u[6]
    ė_θ₃ = θ̇₃d - θ̇₃a

    ė = [ė_θ₁;ė_θ₃]

    τ = p.Kp_f*e + p.Kd_f*ė
end
function Controller_S(u,p,t)
    # Force
    Fd = p.A*sin(p.w*t + p.ϕ)
    Fa = p.kₛ*u[3]
    e_F = Fd - Fa

    # Velecity at Hip
    Vd = p.Vhipd
    Va = Vhip(p,u)
    e_V = Vd - Va

    e = [e_F;e_V]

    # Force derivative
    Ḟd = p.A*p.w*cos(p.w*t + p.ϕ)
    Ḟa = p.kₛ*u[6]
    ė_F = Ḟd - Ḟa

    # Acceteration at CM
    V̇d = 0
    V̇a = Ahip(p,u)
    ė_V = V̇d - V̇a

    ė = [ė_F;ė_V]

    τ = p.Kp_s*e + p.Kd_s*ė
end

##### Flight Functions #####
function FlightTimeApprox(u,p)
    yf = p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2])
    if((u[8]^2-2*p.g*(yf-u[4]))>=0)
        t = (u[8] + sqrt(u[8]^2-2*p.g*(yf-u[4])))/p.g
    else
        t = 0
    end
end
function TrajPlan(u,p)
    tf = p.tf
    M_t =[1   0    0    0;
          0   1    0    0;
          1   tf   tf^2 tf^3;
          0   1    2*tf 3*tf^2]
    M_IC = [u[1:2]'; u[5:6]'; p.setpoint]
     C = M_t\M_IC
end

function Traj(p,t)
     if p.tf < t
          t = p.tf
     end

     C = p.C
     θd = [1 t t^2 t^3;
           0 1 2*t 3*t^2]*C
end
function getFootPos(u,p)
    x_foot = [u[3] + p.Lₐ*sin(u[1]) + p.Lᵦ*sin(u[1]+ u[2]);
              u[4] - p.Lₐ*cos(u[1]) - p.Lᵦ*cos(u[1]+ u[2])]
end
function getFootVel(u,p)
    v_foot = [u[7] + u[5]*p.Lₐ*cos(u[1]) + (u[5]+ u[6])*p.Lᵦ*cos(u[1]+ u[2]);
              u[8] + u[5]*p.Lₐ*sin(u[1]) + (u[5]+ u[6])*p.Lᵦ*sin(u[1]+ u[2])]
end
function getFootAc(u,p)
    v_foot = [u[7] + u[5]*p.Lₐ*cos(u[1]) + (u[5]+ u[6])*p.Lᵦ*cos(u[1]+ u[2]);
              u[8] + u[5]*p.Lₐ*sin(u[1]) + (u[5]+ u[6])*p.Lᵦ*sin(u[1]+ u[2])]
end
##### Stance Functions #####
function Vcm(p,u)
    u = Aug_S(p,u,0)
    J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
         ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
    Pcm = J*u[5:8]
    ẋcm = Pcm[1]
end
function Acm(p,u)
    u = Aug_S(p,u,0)
    J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
         ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
    
    Vⱼ= [-u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])) + -u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*sin(u[1]+u[2]);
          u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])) +  u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*cos(u[1]+u[2])]/(p.m₁ + p.m₂)
    
    Pcm = J*u[5:8] + Vⱼ
    ẋcm = Pcm[1]
end
function Vhip(p,u)
    ẋcm = Aug_S(p,u,0)[7]
end
function Ahip(p,u)
    u = Aug_S(p,u,0)
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

    Phip = -M\(V)
    ẍhip = Phip[3]
end
function Aug_S(p,u,x₀)
    C = [1 0 0;
         0 1 0;
        -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2])    0;
        -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])    1;]
      
    posO = [x₀ + -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])
                  p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1]) + u[3]]

    u = [u[1:2]; posO; C*u[4:6]; u[3]; u[6]]
end

##### State Mapping #####
function q_F(u_all, u)
    u_temp = [reshape(u[1],(1,8)) 0 0]
    for i in 2:length(u)
        u_temp = [u_temp; reshape(u[i],(1,8)) 0 0]
    end
    u_all = [u_all; u_temp]
end
function q_S(u_all, p, u, x₀)
    u_temp = reshape(Aug_S(p,u[1],x₀),(1,10))
    for i in 2:length(u)
        u_temp = [u_temp; reshape(Aug_S(p,u[i],x₀),(1,10))]
    end
    u_all = [u_all; u_temp]
end

### Convert to Cartesian ###
function X_Y(u,p)
    n = length(u[:,1])

    xₒ = u[:,3]
    yₒ = u[:,4]
    u₁ = u[:,1]
    u₂ = u[:,2] + u[:,1]
    yₛ = u[:,9]

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
    cₛ::Float64
    g::Float64

    C::Matrix{Float64}
    tf::Float64
    setpoint::Matrix{Float64}
    Kp_f::Matrix{Float64}
    Kd_f::Matrix{Float64}

    A::Float64
    w::Float64
    ϕ::Float64
    Vhipd::Float64
    Kp_s::Matrix{Float64}
    Kd_s::Matrix{Float64}
end

##################### Start of Code #####################
### Create System ###
p = System(3,.5,.135,.14,.02,.1016,.1,.1,10000,1,9.81,zeros(4, 2),1,zeros(2, 2),zeros(2, 2),zeros(2, 2),0,0,0,0,zeros(2, 2),zeros(2, 2))

## Start in Flight ##
x₀ = [.9 -1.6 0 .4 -.3 .8 -.5 0]
x₀ = [.6 -1.2 .24 .34 0 0 0 0]

### Set up Fight Controllers ###
p.setpoint = [0.6 -1.1; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### Set up Controller ###
p.A = 20
p.w = 50*2*pi
p.Vhipd = 1.5
p.Kp_s = [.185 .21; .2 .1]
p.Kd_s = [.0 0; 0 0]

####### ODE Solver #######
tEnd = 3
affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3] + (u[6]<0)*-100
c_impact(u,t,integrator) = getFootPos(u,integrator.p)[2] + (getFootVel(u,integrator.p)[2]>=0)*+100
cb_lift = ContinuousCallback(c_lift,affect!)
cb_impact = ContinuousCallback(c_impact,affect!)

### Initial Conditions ###
t_cur = 0
t_all = [t_cur]
u_all = [x₀ 0 0]
sw = [1]

###### Run Simulation ######
while t_cur < tEnd -.01
    global t_cur, x₀, u_all, t_all, sw
    tspan = (0, tEnd - t_cur)

    p.tf = FlightTimeApprox(x₀,p);
    p.C = TrajPlan(x₀,p);
    prob = ODEProblem(ode_F_fn!,x₀,tspan,p);
    sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = .01);

    if length(sol.u[:,1]) >= 2 && sol.u[end-1,:] != sol.u[end,:]
        u_all = q_F(u_all, sol.u[2:end,:])
        t_all = [t_all; ones(length(sol.t)-1,1)*t_cur + sol.t[2:end,:]]
        t_cur = t_all[end]
        sw = [sw; length(u_all[:,1])]
    elseif length(sol.u[:,1]) >= 3
        u_all = q_F(u_all, sol.u[2:end-1,:])
        t_all = [t_all; ones(length(sol.t)-2,1)*t_cur + sol.t[2:end-1,:]]
        t_cur = t_all[end]
        sw = [sw; length(u_all[:,1])]
    else
        sw = [sw; length(u_all[:,1])]
    end

## Rosenbrock23()
## Rosas5P(): increase effecency; slower computation
    if t_cur < tEnd -.01
        x₀ = [u_all[end,1:2]; getFootPos(u_all[end,:],p)[2]; u_all[end,5:6]; getFootVel(u_all[end,:],p)[2]]
        x_foot = getFootPos(u_all[end,:],p)[1]

        tspan = (0, tEnd - t_cur)
        prob = ODEProblem(ode_S_fn!,x₀,tspan,p)
        sol = solve(prob,AutoTsit5(Rodas4P()),callback=cb_lift, abstol=1e-9,reltol=1e-9, saveat = .01)

        if length(sol.u[:,1]) >= 2 && sol.u[end-1,:] != sol.u[end,:]
            u_all = q_S(u_all, p, sol.u[2:end,:], x_foot)
            t_all = [t_all; ones(length(sol.t)-1,1)*t_cur + sol.t[2:end,:]]
            t_cur = t_all[end]
            sw = [sw; length(u_all[:,1])]
        elseif length(sol.u[:,1]) >= 3
            u_all = q_S(u_all, p, sol.u[2:end-1,:], x_foot)
            t_all = [t_all; ones(length(sol.t)-2,1)*t_cur + sol.t[2:end-1,:]]
            t_cur = t_all[end]
            sw = [sw; length(u_all[:,1])]
        else
            sw = [sw; length(u_all[:,1])]
        end
        x₀ = u_all[end,1:8]
    end
end

##### Export Data #####
df1 = DataFrame(t = t_all[:,1], θ₁=u_all[:,1], θ₃=u_all[:,2], x=u_all[:,3], y=u_all[:,4], θ̇₁=u_all[:,5], θ̇₃=u_all[:,6], ẋ=u_all[:,7], ẏ=u_all[:,8], yₛ=u_all[:,9], ẏₛ =u_all[:,10])
df2 = DataFrame(tf = [p.tf, missing], Kp_f1 = p.Kp_f[:,1], Kp_f2 = p.Kp_f[:,2], Kd_f1 = p.Kd_f[:,1], Kd_f2 = p.Kd_f[:,2], A=[p.A, missing], w=[p.w, missing], ϕ=[p.ϕ, missing], Vhipd=[p.Vhipd, missing], Kp_s1 = p.Kp_s[:,1], Kp_s2 = p.Kp_s[:,2], Kd_s1 = p.Kd_s[:,1], Kd_s2 = p.Kd_s[:,2])
df3 = DataFrame(m₁=p.m₁,m₂=p.m₂,Lₐ=p.Lₐ,Lᵦ=p.Lᵦ,L₁=p.L₁,L₂=p.L₂,c₁=p.c₁,c₂=p.c₂,kₛ=p.kₛ,cₛ=p.cₛ,g=p.g)

cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData")
a = readdir()
num = length(a)+1

XLSX.writetable("z_DataSet" * string(num)*".xlsx", overwrite=true, 
    OUTPUT_DATA=(collect(DataFrames.eachcol(df1)), DataFrames.names(df1)),
    CONTROL_PARAM=(collect(DataFrames.eachcol(df2)), DataFrames.names(df2)),
    SYSTEM_PARAM=(collect(DataFrames.eachcol(df3)), DataFrames.names(df3)), 
)
#=
open("test1.txt", "a") do file
    write(file, string(y))
end
=#
