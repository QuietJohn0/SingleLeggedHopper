################## Double Pendulum XLSX Plotter 1.0 ##################

################### Imported Librarys ###################
#using DifferentialEquations
#using LinearAlgebra
using StaticArrays
using XLSX
using DataFrames
using Plots
#using PlotlyJS


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

### Fin Functions ###
function Findsw(u_all)
    sw = [1]
    for i in 2:length(u_all[:,1])
        if u_all[i,9] != u_all[i-1,9] && u_all[i,9]*u_all[i-1,9] == 0
            sw = [sw i]
        end
    end
    sw = [sw length(u_all[:,1])+1]
end

### Plotting ###
function PlotTransitionLines(t_all,u, title,sw)
    plot(t_all, u, title=title, legend=false)
    for i in 1:length(sw)-1
        j = sw[i]
        k = sw[i+1]
        if i%2 == 0
            plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
            plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
        else
            plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)
            plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)        
        end
    end
    plot!([t_all[end] t_all[end]], [findmin(u)[1] findmax(u)[1]], lc=:green, lw=1)
end
function PlotTransitionLines!(t_all,u, title,sw)
    plot!(t_all, u, title=title, legend=false)
    for i in 1:length(sw)-1
        j = sw[i]
        k = sw[i+1]
        if i%2 == 0
            plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
            plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
        else
            plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)
            plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)        
        end
    end
    plot!([t_all[end] t_all[end]], [findmin(u)[1] findmax(u)[1]], lc=:green, lw=1)
end
function PlotTransitionPoints(ux,uy, title,sw)
    plot(ux, uy, title=title, legend=false)
    for i in 1:length(sw)-1
        j = sw[i]
        k = sw[i+1]
        if i%2 == 0

            plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :red, shape = :circle, markersize = 2)
            plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :red, shape = :circle, markersize = 2)
        else
            #Flight
            plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :green, shape = :circle, markersize = 2)
            plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :green, shape = :circle, markersize = 2)        
        end
    end
    plot!([ux[end] ux[end]], [uy[end] uy[end]], markercolors = :black, shape = :circle, markersize = 2)
end
function PlotTransitionPoints!(ux,uy, title,sw)
    plot!(ux, uy, title=title, legend=false)
    for i in 1:length(sw)-1
        j = sw[i]
        k = sw[i+1]
        if i%2 == 0
            plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :red, shape = :circle, markersize = 2)
            plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :red, shape = :circle, markersize = 2)
        else
            plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :green, shape = :circle, markersize = 2)
            plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :green, shape = :circle, markersize = 2)        
        end
    end
    plot!([ux[end] ux[end]], [uy[end] uy[end]], markercolors = :black, shape = :circle, markersize = 2)
end
function PlotSegmentCompare(t_all,u,y, title,sw_1, sw_2)
    t = t_all[sw_1:sw_2]
    u = u[sw_1:sw_2]
    plot(t, u, title=title, legend=false, lc=blue, lw=.7)
    plot!(t, y, lc=black, ls=:dash, lw=.7)
end

### Control Error Calc and Plots###
function PlotStanceGRF(t_all,u,p, title,sw)
    PlotTransitionLines(t_all, -p.kₛ*u_all[:,9] - p.cₛ*u[:,10], title, sw)
    Fd = p.A*sin.(p.w.*t_all+ones(size(t_all))*p.ϕ)
    PlotTransitionLines!(t_all, Fd, title, sw)
end
function PlotStanceGRFDev(t_all,u,p, title,sw)
    PlotTransitionLines(t_all, -p.kₛ*u[:,10], title, sw)
    Ḟd = p.A*p.w.*cos.(p.w.*t_all+ones(size(t_all))*p.ϕ)
    PlotTransitionLines!(t_all, Ḟd, title, sw)
end
function PlotStanceV(t_all,u,p, title,sw)
    PlotTransitionLines(t_all, u_all[:,7], title, sw)
    plot!([t_all[1] t_all[end]], p.Vhipd*ones(1,2))
end
function PlotStanceA(t_all,u_all,p, title,sw)
    display(u_all)
    PlotTransitionLines(t_all, Ahipl(p,u_all), title, sw)
    plot!([t_all[1] t_all[end]], 0*ones(1,2))
end


### Torque ###
function findTorques(u_all,p,t_all, sw)
    τ = zeros(2,0)
    for i in 1:length(sw)-1
        j = sw[i]
        k = sw[i+1]-1
        if i%2 == 0
            #Stance
            for l in j:k
                τ = [τ Controller_S(u_all[l,:],p,t_all[l]-t_all[j])]
            end
        else
            #Flight
            print(t_all[j])
            #p.C = TrajPlan(u_all[j,1:8],p.tf, p.setpoint)
            for l in j:k
                τ = [τ Controller_F(u_all[l,:],p,t_all[l]-t_all[j])]
            end
        end
    end
    τ
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

    # Acceteration at CM, θ̇₃
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

    # Velecity at CM
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

##### Stance Functions #####
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
function Ahipl(p,u_all)
    ẍhip = Ahip(p,Aug(p,u_all[1,:],0))
    for i in 2:length(u_all[:,end])
        ẍhip = [ẍhip; Ahip(p,Aug(p,u_all[i,:],0))]
    end
end

##### Flight Functions #####
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
    M_IC = [u[1:2]; u[5:6]; setpoint]
    C = M_t\M_IC
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


##### Export Data #####
cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData")
a = readdir()

# Open the Excel file
xf = XLSX.readxlsx(a[end])
sheetNames = XLSX.sheetnames(xf)

# Extract from sheet 1
sh = xf[sheetNames[1]]
t_all = sh["A"][2:end]
u_all = sh["B:K"][2:end,:]

# Extract from sheet 2
sh = xf[sheetNames[2]]
p2 = sh["A:M"][2:end,:]

# Extract from sheet 3
sh = xf[sheetNames[3]]
p1 = sh["A:K"][2,:]

# Build System
p = System(p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8],p1[9],p1[10],p1[11],
        zeros(4, 2),p2[1,1],zeros(2, 2),[p2[:,2] p2[:,3]],[p2[:,4] p2[:,5]],p2[1,6],p2[1,7],p2[1,8],p2[1,9],[p2[:,10] p2[:,11]],[p2[:,12] p2[:,13]])


# Generate some data
x,y = X_Y(u_all,p)
sw = Findsw(u_all)
τ = findTorques(u_all,p,t_all,sw)
plotlyjs()

# Create a grid of subplots
p1 = PlotTransitionLines(t_all, u_all[:,1], "θ₁", sw)
p2 = PlotTransitionLines(t_all, u_all[:,2], "θ₃", sw)
p3 = PlotTransitionLines(t_all, u_all[:,3], "x", sw)
p4 = PlotTransitionLines(t_all, u_all[:,4], "y", sw)
p5 = PlotTransitionLines(t_all, u_all[:,5], "θ̇₁", sw)
p6 = PlotTransitionLines(t_all, u_all[:,6], "θ̇₃", sw)
p7 = PlotTransitionLines(t_all, u_all[:,7], "ẋ", sw)
p8 = PlotTransitionLines(t_all, u_all[:,8], "ẏ", sw)
sp1 = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(2,4))
display(sp1)

p1 = PlotTransitionPoints(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)
p2 = PlotTransitionPoints(u_all[:,2],u_all[:,6], "θ₃ v.s. θ̇₃",sw)
p3 = PlotTransitionPoints(u_all[:,3],u_all[:,4], "x v.s. y",sw)
p4 = PlotTransitionLines(t_all, p.kₛ*u_all[:,9], "grf", sw)


p5 = PlotTransitionPoints(u_all[:,1], u_all[:,1] + u_all[:,2], "θ₁ v.s. θ₁+θ₃", sw)
p6 = PlotTransitionPoints(t_all, τ[1,:], "Torque, τ", sw)
PlotTransitionPoints!(t_all, τ[2,:], "Torque, τ", sw)
sp2 = plot(p1, p2, p3, p4, p5, p6, layout=(2,3))
display(sp2)

PlotTransitionPoints(t_all, u_all[:,8], "Vhipd VS Vhipd", sw)
plot!([t_all[1]; t_all[end]],[p.Vhipd; p.Vhipd],lc=:black, ls=:dash,lw = 2)

### Stance Controls ###
display(PlotStanceGRF(t_all,u_all,p, "GRF",sw))
display(PlotStanceGRFDev(t_all,u_all,p, "GṘF",sw))
display(PlotStanceV(t_all,u_all,p, "Vhid",sw))
display(PlotStanceA(t_all,u_all,p, "Ahip",sw))


#p1 = PlotSegmentCompare(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)

τ = findTorques(u_all,p,t_all,sw)

#=############## Phase Plot ###############
##### X V.S Y #####
x,y = X_Y(u_all,p)
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
display(plt)=#
