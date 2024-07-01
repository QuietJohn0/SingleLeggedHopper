################## Double Pendulum 3.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays
using XLSX
using DataFrames

####################### Functions #######################
### ODE Solver ###
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

function Controller_S(u,p,t)
    # Force
    Fd = p.A*sin(p.w*t + p.ϕ)
    Fa = -p.kₛ*u[3] - p.cₛ*u[6]
    e_F = Fd - Fa

    # Velecity at Hip
    Vd = p.Vhipd
    Va = Vhip(p,u)
    e_V = Vd - Va

    e = [e_F;e_V]

    # Force derivative
    Ḟd = p.A*p.w*cos(p.w*t + p.ϕ)
    Ḟa = -p.kₛ*u[6] - p.cₛ*u[6]
    ė_F = Ḟd - Ḟa

    # Acceteration at CM
    V̇d = 0
    V̇a = Ahip(p,u)
    ė_V = V̇d - V̇a

    ė = [ė_F;ė_V]

    τ = p.Kp_s*e + p.Kd_s*ė
end

function Vcm(p,u)
     u = Aug(p,u,0)
     J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
          ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
     Pcm = J*u[5:8]
     ẋcm = Pcm[1]
end
function Acm(p,u)

     u = Aug(p,u,0)
     J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₁ + p.m₂)    1   0;
          ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₁ + p.m₂)    0   1]
     
     Vⱼ= [-u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])) + -u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*sin(u[1]+u[2]);
           u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])) +  u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*cos(u[1]+u[2])]/(p.m₁ + p.m₂)
     
     Pcm = J*u[5:8] + Vⱼ
     ẋcm = Pcm[1]
end
function Vhip(p,u)
    ẋcm = Aug(p,u,0)[7]
end
function Ahip(p,u)
    u = Aug(p,u,0)
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
function Ahip!(p,u_all)
    ẍhip = Ahip(p,Aug(p,u_all[1,:],0))
    for i in 2:length(u_all[:,end])
        ẍhip = [ẍhip; Ahip(p,Aug(p,u_all[i,:],0))]
    end
end

function getFootVel(u,p)
    v_foot = [u[7] + u[5]*p.Lₐ*cos(u[1]) + (u[5]+ u[6])*p.Lᵦ*cos(u[1]+ u[2]);
              u[8] + u[5]*p.Lₐ*sin(u[1]) + (u[5]+ u[6])*p.Lᵦ*sin(u[1]+ u[2])]
end
function PhiShiftCalc(x₀,p)
    if -(p.kₛ*x₀[3] + p.cₛ*x₀[6])>p.A
        p.A
    elseif -(p.kₛ*x₀[3] + p.cₛ*x₀[6])<0
        0
    else
        asin(-(p.kₛ*x₀[3] + p.cₛ*x₀[6])/p.A)
    end
end

##### State Mapping #####
function q(u_all, p, u, x₀)
    u_temp = reshape(Aug(p,u[1],x₀),(1,10))
    for i in 2:length(u)
        u_temp = [u_temp; reshape(Aug(p,u[i],x₀),(1,10))]
    end
    u_all = [u_all; u_temp]
end
function Aug(p,u,x₀)
    if (length(u) == 6)
        out = Aug_S(p,u,x₀)
    elseif length(u) == 8
        out = [u 0 0]
    elseif length(u) == 10
        out = u
    end
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
p = System(3,.5,.135,.14,.02,.1016,.1,.1,20000,200,9.81,zeros(4, 2),1,zeros(2, 2),zeros(2, 2),zeros(2, 2),0,0,0,0,zeros(2, 2),zeros(2, 2))

### Set up Stance Controller ###
p.A = 40
p.w = 20*pi
p.ϕ = 0
p.Vhipd = -.25
p.Kp_s = [1.8 0; 2 0]
p.Kd_s = [.1 0; .1 0]

### Initial Conditions ###
x₀ = [1, -2, 0, 0, 0, 0];
x₀ = [0.6, -1.2,  0.0,  1.29658939587674e-8, -2.583073140309899e-8,  0.0];
x₀ = [0.5988020714847644 -1.098852191182763  0 -3.9740706086530313 5.929846985702695 -1.7834013756061506];
x₀ = [0.6323468141422804 -1.1500911731708694 0.01635024771144429 -3.0252402252359136 4.757556756433081 -1.6350247711444286]

####### ODE Solver #######
tspan = (0, 5);

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3] + (u[6]<0)*-100
cb_lift = ContinuousCallback(c_lift,affect!)

p.ϕ = PhiShiftCalc(x₀,p)
prob = ODEProblem(ode_S_fn!,x₀,tspan,p)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = .001)

##### Interpret Data #####
t_all = [0]


u_all = q(missing, p, sol.u[1:end,:], 0)[2:end,1:end]
t_all = [[]; sol.t[1:end,:]]

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
 end every 10=#
