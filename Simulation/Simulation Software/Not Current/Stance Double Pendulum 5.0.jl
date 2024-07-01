################## Double Pendulum 5.0 ##################

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
     τ =  [Controller_S(u,p,t); 0];

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

### Controller ###
function Controller_S(u,p,t)
    # Feed Forward Torque
    Fd = [p.Fxd; p.A*sin(p.w*t + p.ϕ)]
    τff = FeedForwardTorque(p,u,Fd)

    # Hip y Position
    yd = p.yhipd
    ya = Phip(p,u)[2]
    e_y = yd - ya

    # Hip y Velecity
    Vyd = p.Vyhipd
    Vya = Vhip(p,u)[2]
    e_Vy = Vyd - Vya

    # Hip x Velecity
    Vxd = p.Vxhipd
    Vxa = Vhip(p,u)[1]
    e_Vx = Vxd - Vxa

    # Compute error Matricies
    e = e_y
    ė = [e_Vx; e_Vy]

    τ = τff + p.Kp_s*e + p.Kd_s*ė

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
function Phip(p,u)
    xhip = Aug(p,u,0)[3:4]
end
function Vhip(p,u)
   ẋhip = Aug(p,u,0)[7:8]
end
function APhiGRFCalc(x₀,p)
    u = Aug(p,x₀,0)

    A = (((p.m₁+p.m₂)*p.w*(p.Vyhipd - u[8]))^2 + (p.cₛ*u[10])^2)/((p.m₁+p.m₂)*p.w*(p.Vyhipd - u[8]))
    ϕ = asin(-p.cₛ*u[10]/A)

    (A, ϕ)
end
function FeedForwardTorque(p,u,F)
    J = [p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1])  p.Lᵦ*cos(u[1]+u[2]);
         p.Lᵦ*sin(u[1]+u[2]) + p.Lₐ*sin(u[1])  p.Lᵦ*sin(u[1]+u[2])]

    τ_mass = [p.m₁*p.g*p.L₁*sin(u[1]) + p.m₂*p.g*(p.L₂*sin(u[1]+u[2]) + p.Lₐ*sin(u[1]));
              p.m₂*p.g*p.L₂*sin(u[1]+u[2])]

    τff = -J'*F - τ_mass
end
function VhipCalc(p)
    ## xhop is an approximate distance
    yf = p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2])
    Vyhipd = sqrt(2*p.g*(p.y_hop-p.yhipd))
    Vxhipd = p.dx_hop*p.g/(Vyhipd + sqrt(Vyhipd^2 + 2*p.g*(p.yhipd - yf)))
    (Vxhipd,Vyhipd)
end



##### Transition Functions #####
function getFootPos(u_all,p)
    x_foot = [u_all[3] + p.Lₐ*sin(u_all[1]) + p.Lᵦ*sin(u_all[1]+ u_all[2]);
              u_all[4] - p.Lₐ*cos(u_all[1]) - p.Lᵦ*cos(u_all[1]+ u_all[2])]
end
function getFootVel(u_all,p)
    v_foot = [u_all[7] + u_all[5]*p.Lₐ*cos(u_all[1]) + (u_all[5]+ u_all[6])*p.Lᵦ*cos(u_all[1]+ u_all[2]);
              u_all[8] + u_all[5]*p.Lₐ*sin(u_all[1]) + (u_all[5]+ u_all[6])*p.Lᵦ*sin(u_all[1]+ u_all[2])]
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

##### State Mapping #####
    function q(u_all, p, u, x_foot)
        u_temp = reshape(Aug(p,u[1],x_foot),(1,10))
        for i in 2:length(u)
            u_temp = [u_temp; reshape(Aug(p,u[i],x_foot),(1,10))]
        end
        u_all = [u_all; u_temp]
    end
    function Aug(p,u,x_foot)
        if (length(u) == 6)
            out = Aug_S(p,u,x_foot)
        elseif length(u) == 8
            out = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] getFootPos(u,p)[2] getFootVel(u,p)[2]]
        elseif length(u) == 10
            out = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] u[9] u[10]]
        end
    end
    function invAug(p,u,n)
        if (n == 6)
            x₀ = [u[1] u[2] getFootPos(u,p)[2] u[5] u[6] getFootVel(u,p)[2]]
        elseif n == 8
            x₀ = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8]]
        elseif n == 10
            x₀ = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] getFootPos(u,p)[2] getFootVel(u,p)[2]]
        end
        x_foot = getFootPos(u,p)[1]

        (x₀, x_foot)
    end
    function Aug_S(p,u,x_foot)
        C = [1 0 0;
            0 1 0;
            -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2])    0;
            -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])    1;]
        posȮ = C*u[4:6]
        
        posO = [x_foot + -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])
                        p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1]) + u[3]]

        u = [u[1] u[2] posO[1] posO[2] posȮ[1] posȮ[2] posȮ[3] posȮ[4] u[3] u[6]]
    end
### END State Mapping ###

##### Importing and Exporting Data #####
    function ExportData(t_all,u_all,p,name)
        df1 = DataFrame(t = t_all[:,1], θ₁=u_all[:,1], θ₃=u_all[:,2], x=u_all[:,3], y=u_all[:,4], θ̇₁=u_all[:,5], θ̇₃=u_all[:,6], ẋ=u_all[:,7], ẏ=u_all[:,8], yₛ=u_all[:,9], ẏₛ =u_all[:,10])
        df2 = DataFrame(setpoint1 = p.setpoint[:,1], setpoint2 = p.setpoint[:,2] , Kp_f1 = p.Kp_f[:,1], Kp_f2 = p.Kp_f[:,2], Kd_f1 = p.Kd_f[:,1], Kd_f2 = p.Kd_f[:,2], 
                        w=[p.w, missing], Fxd=[p.Fxd, missing], yhipd=[p.yhipd, missing], Kp_s = p.Kp_s[:,1], Kd_s1 = p.Kd_s[:,1], Kd_s2 = p.Kd_s[:,2],
                        dx_hop=[p.dx_hop, missing], y_hop=[p.y_hop, missing])
        df3 = DataFrame(m₁=p.m₁,m₂=p.m₂,Lₐ=p.Lₐ,Lᵦ=p.Lᵦ,L₁=p.L₁,L₂=p.L₂,c₁=p.c₁,c₂=p.c₂,kₛ=p.kₛ,cₛ=p.cₛ,g=p.g)

        cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData")
        a = readdir()
        num = length(a)+1

        XLSX.writetable(name * string(num)*".xlsx", overwrite=true, 
            OUTPUT_DATA=(collect(DataFrames.eachcol(df1)), DataFrames.names(df1)),
            CONTROL_PARAM=(collect(DataFrames.eachcol(df2)), DataFrames.names(df2)),
            SYSTEM_PARAM=(collect(DataFrames.eachcol(df3)), DataFrames.names(df3)), 
        )
    end
    function ImportData()
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
        p2 = sh["A:N"][2:end,:]

        # Extract from sheet 3
        sh = xf[sheetNames[3]]
        p1 = sh["A:K"][2,:]

        # Build System
        p = System(p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8],p1[9],p1[10],p1[11],
                zeros(4, 2),0,[p2[:,1] p2[:,2]],[p2[:,3] p2[:,4]],[p2[:,5] p2[:,6]],
                0,p2[1,7],0,p2[1,8],p2[1,9],0,0,p2[:,10],[p2[:,11] p2[:,12]],
                p2[1,13],p2[1,14])

        (t_all,u_all,p,a[end])
    end
### END Importing and Exporting Data ###

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
    Fxd::Float64
    yhipd::Float64
    Vxhipd::Float64
    Vyhipd::Float64
    Kp_s::Vector{Float64}
    Kd_s::Matrix{Float64}

    dx_hop::Float64
    y_hop::Float64
end

##################### Start of Code #####################
### Create System ###
p = System(1.6,.2,.135,.164,.013,.085,.1,2,20000,20,9.81,
            zeros(4, 2),0,zeros(2, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,0,[0; 0],zeros(2, 2),
            0,0)

### Set up Fight Controller ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### Set up Stance Controller ###
p.w = 50*pi
p.Fxd = 0
p.yhipd = (p.Lᵦ + p.Lₐ + p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2]))/2
p.Kp_s = [200; 70]*0
p.Kd_s = [0 0; 40 0]

#p.Kp_s = [2000; 700]*0
#p.Kd_s = [4 100; 200 25]

### System Goals ###
p.dx_hop = -.3
p.y_hop = .32
(p.Vxhipd,p.Vyhipd) = VhipCalc(p)

### Initial Conditions ###
u_all = [0.60	-1.03	0.20	0.26	-3.75	4.83	-0.23	-1.20	0.00	-1.56]
(u_all,x_foot) = invAug(p,u_all,10)

####### ODE Solver #######
tspan = (0, 5);

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3] + (u[6]<0)*-100 #-integrator.p.kₛ*u[3] - integrator.p.cₛ*u[6] + (u[6]<0)*100
cb_lift = ContinuousCallback(c_lift,affect!)

(x₀,x_foot) = invAug(p,u_all[end,:],6)

(p.A,p.ϕ) = APhiGRFCalc(x₀,p)

prob = ODEProblem(ode_S_fn!,x₀,tspan,p)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = .001)

##### Interpret Data #####
u_all = q(u_all, p, sol.u[2:end,:], x_foot)
t_all = [[0]; sol.t[2:end,:]]

##### Export Data #####
name = "z_DataSet"
ExportData(t_all,u_all,p,name)

