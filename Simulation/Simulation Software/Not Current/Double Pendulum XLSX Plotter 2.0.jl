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

### Find Switches Functions ###
function Findsw(u_all)
    sw = [0]
    if (u_all[1,9] < 0)
        sw = [sw 1]
    end

    for i in 2:length(u_all[:,1])
        if (u_all[i,9] < 0 && u_all[i-1,9] >= 0)
            sw = [sw i]
        elseif (u_all[i,9] > 0 && u_all[i-1,9] <= 0)
            sw = [sw i-1]
        end
    end
    sw = [sw length(u_all[:,1])]
end

##### State Transition Plots #####
    function PlotTransitionLines(t_all,u, title,sw)
        plot(t_all, u, title=title, legend=false)
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
            elseif j != k
                plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)        
            end
        end
        plot!([t_all[end] t_all[end]], [findmin(u)[1] findmax(u)[1]], lc=:green, lw=1)
    end
    function PlotTransitionLines!(t_all,u, title,sw)
        plot!(t_all, u, title=title, legend=false)
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:red, ls=:dash, lw=.7)
            elseif j != k
                plot!([t_all[j], t_all[j]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [findmin(u)[1], findmax(u)[1]], lc=:green, ls=:dash, lw=.7)        
            end
        end
        plot!([t_all[end] t_all[end]], [findmin(u)[1] findmax(u)[1]], lc=:green, lw=1)
    end
    function PlotTransitionPoints(ux,uy, title,sw)
        plot(ux, uy, title=title, legend=false)
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0

                plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :red, shape = :circle, markersize = 2)
                plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :red, shape = :circle, markersize = 2)
            elseif j != k
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
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                plot!([ux[j], ux[j]], [uy[j], uy[j]], markercolors = :red, shape = :circle, markersize = 2)
                plot!([ux[k-1], ux[k-1]], [uy[k-1], uy[k-1]], markercolors = :red, shape = :circle, markersize = 2)
            elseif j != k
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
### END State Transition Plots ###

##### Stance Controller Plots #####
    function PlotStanceGRF(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, -p.kₛ*u_all[:,9] - p.cₛ*u[:,10], title, sw)
        Fd = p.A*sin.(p.w.*t_all+ones(size(t_all))*p.ϕ)
        PlotTransitionLines!(t_all, Fd, title, sw)
    end
    function PlotStancePy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,4], title, sw)
        plot!([t_all[1], t_all[end]], p.yhipd*[1,1])
    end
    function PlotStanceVx(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,7], title, sw)
        plot!([t_all[1], t_all[end]], p.Vhipd*[1,1])
    end
    function PlotStanceVy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,8], title, sw)
        plot!([t_all[1], t_all[end]], [0,0])
    end
### END Stance Controller Plots ###




### Torque ###
function findTorques(u_all,p,t_all, sw)
    τ = zeros(2,0)
    for i in 1:length(sw)-1
        j = sw[i]+1
        k = sw[i+1]
        if i%2 == 0
            #Stance
            p.ϕ = PhiShiftCalc(u_all[j,:],p)
            for l in j:k
                τ = [τ Controller_S(u_all[l,:],p,t_all[l]-t_all[j])]
            end
        else
            #Flight
            p.tf = FlightTimeApprox(u_all[j,:],p);
            for l in j:k
                p.C = TrajPlan(u_all[j,:],p);
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

        #  Position Derivative, θ̇₃
        θ̇₃d = θd[2,2]
        θ̇₃a = u[6]
        ė_θ₃ = θ̇₃d - θ̇₃a

        ė = [ė_θ₁;ė_θ₃]

        τ = p.Kp_f*e + p.Kd_f*ė
    end
    function Controller_S(u,p,t)
        # Feed Forward Torque
        Fd = [p.Fxd; p.A*sin(p.w*t + p.ϕ)]
        τff = FeedForwardTorque(p,u,Fd)

        # Hip y Position
        yd = p.yhipd
        ya = Phip(p,u)[2]
        e_y = yd - ya

        # Hip y Velecity
        Vyd = 0
        Vya = Vhip(p,u)[2]
        e_Vy = Vyd - Vya

        # Hip x Velecity
        Vxd = p.Vhipd
        Vxa = Vhip(p,u)[1]
        e_Vx = Vxd - Vxa

        # Compute error Matricies
        e = e_y
        ė = [e_Vx; e_Vy]

        τ = τff + p.Kp_s*e + p.Kd_s*ė

    end
###### END Controllers ######

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
    function PhiShiftCalc(x₀,p)
        u = Aug(p,x₀,0)
        if -(p.kₛ*u[9] + p.cₛ*u[10])>p.A
            0
        elseif -(p.kₛ*u[9] + p.cₛ*u[10])<0
            0
        else
            asin(-(p.kₛ*u[9] + p.cₛ*u[10])/p.A)
        end
    end
    function FeedForwardTorque(p,u,F)
        J = [p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1])  p.Lᵦ*cos(u[1]+u[2]);
             p.Lᵦ*sin(u[1]+u[2]) + p.Lₐ*sin(u[1])  p.Lᵦ*sin(u[1]+u[2])]
    
        τ_mass = [p.m₁*p.g*p.L₁*sin(u[1]) + p.m₂*p.g*(p.L₂*sin(u[1]+u[2]) + p.Lₐ*sin(u[1]));
                  p.m₂*p.g*p.L₂*sin(u[1]+u[2])]
    
        τff = -J'*F - τ_mass
    end
### END Stance Functions ###

##### Flight Functions #####
    function FlightTimeApprox(u,p)
        yf = p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2])
        if((u[8]^2-2*p.g*(yf-u[4]))>=0)
            t = (u[8] + sqrt(u[8]^2-2*p.g*(yf-u[4])))/p.g
        else
            t = 0.005
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
### END Flight Functions ###

##### Transition Functions #####
    function getFootPos(u_all,p)
        x_foot = [u_all[3] + p.Lₐ*sin(u_all[1]) + p.Lᵦ*sin(u_all[1]+ u_all[2]);
                u_all[4] - p.Lₐ*cos(u_all[1]) - p.Lᵦ*cos(u_all[1]+ u_all[2])]
    end
    function getFootVel(u_all,p)
        v_foot = [u_all[7] + u_all[5]*p.Lₐ*cos(u_all[1]) + (u_all[5]+ u_all[6])*p.Lᵦ*cos(u_all[1]+ u_all[2]);
                u_all[8] + u_all[5]*p.Lₐ*sin(u_all[1]) + (u_all[5]+ u_all[6])*p.Lᵦ*sin(u_all[1]+ u_all[2])]
    end
### END Transition Functions ###

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
        df2 = DataFrame(setpoint1 = p.setpoint[:,1], setpoint2 = p.setpoint[:,2] , Kp_f1 = p.Kp_f[:,1], Kp_f2 = p.Kp_f[:,2], Kd_f1 = p.Kd_f[:,1], Kd_f2 = p.Kd_f[:,2], A=[p.A, missing], w=[p.w, missing], ϕ=[p.ϕ, missing], Fxd=[p.Fxd, missing], yhipd=[p.yhipd, missing], Vhipd=[p.Vhipd, missing], Kp_s = p.Kp_s[:,1], Kd_s1 = p.Kd_s[:,1], Kd_s2 = p.Kd_s[:,2])
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
        t_all = Float64.(sh["A"][2:end])
        u_all = Float64.(sh["B:K"][2:end,:])

        # Extract from sheet 2
        sh = xf[sheetNames[2]]
        p2 = sh["A:O"][2:end,:]

        # Extract from sheet 3
        sh = xf[sheetNames[3]]
        p1 = sh["A:K"][2,:]

        # Build System
        p = System(p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8],p1[9],p1[10],p1[11],
                zeros(4, 2),0,[p2[:,1] p2[:,2]],[p2[:,3] p2[:,4]],[p2[:,5] p2[:,6]],p2[1,7],p2[1,8],p2[1,9],p2[1,10],p2[1,11],p2[1,12],p2[:,13],[p2[:,14] p2[:,15]])

        (t_all,u_all,p)
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
    Vhipd::Float64
    Kp_s::Vector{Float64}
    Kd_s::Matrix{Float64}
end

##### Import Data #####
(t_all,u_all,p) = ImportData()

# Generate some data
x,y = X_Y(u_all,p)
sw = Findsw(u_all)
τ = findTorques(u_all,p,t_all,sw)
plotlyjs()

##### Subplot 1: State Space #####
p1 = PlotTransitionLines(t_all, u_all[:,1], "θ₁", sw)
p2 = PlotTransitionLines(t_all, u_all[:,2], "θ₃", sw)
p3 = PlotTransitionLines(t_all, u_all[:,3], "x", sw)
p4 = PlotTransitionLines(t_all, u_all[:,4], "y", sw)
p5 = PlotTransitionLines(t_all, u_all[:,5], "θ̇₁", sw)
p6 = PlotTransitionLines(t_all, u_all[:,6], "θ̇₃", sw)
p7 = PlotTransitionLines(t_all, u_all[:,7], "ẋ", sw)
p8 = PlotTransitionLines(t_all, u_all[:,8], "ẏ", sw)
sp1 = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(2,4))
#display(sp1)

##### Subplot 2: TBD #####
p1 = PlotTransitionPoints(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)
p2 = PlotTransitionPoints(u_all[:,2],u_all[:,6], "θ₃ v.s. θ̇₃",sw)
p3 = PlotTransitionPoints(u_all[:,3],u_all[:,4], "x v.s. y",sw)
p4 = PlotTransitionLines(t_all, p.kₛ*u_all[:,9], "grf", sw)
p5 = PlotTransitionPoints(u_all[:,1], u_all[:,1] + u_all[:,2], "θ₁ v.s. θ₁+θ₃", sw)
p6 = PlotTransitionPoints(t_all, τ[1,:], "Torque, τ", sw)
p6 = PlotTransitionPoints!(t_all, τ[2,:], "Torque, τ", sw)
sp2 = plot(p1, p2, p3, p4, p5, p6, layout=(2,3))
#display(sp2)

##### Subplot 3: Stance Controls #####
s1 = PlotStanceGRF(t_all,u_all,p, "GRF",sw)
s2 = PlotStancePy(t_all,u_all,p, "Pyhid",sw)
s3 = PlotStanceVx(t_all,u_all,p, "Vxhid",sw)
s4 = PlotStanceVy(t_all,u_all,p, "Vyhid",sw)
sp3 = plot(s1, s2, s3, s4, layout=(2,2))
display(sp3)

#p1 = PlotSegmentCompare(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)
