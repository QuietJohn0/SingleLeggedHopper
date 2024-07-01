################## Double Pendulum XLSX Animator 1.0 ##################

################### Imported Librarys ###################
using StaticArrays
using XLSX
using DataFrames
using Plots



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

##### Animating #####
    function PlotState!(x,y, title)
        plot!(x,y,ylim = [-.25,.85],xlim=[-1,1],
            xlabel="X Position (m)", 
            ylabel="Y Position (m)",
            title=title,
            markercolors = [:blue, :red, :blue, :red, :black], 
            markershape = [:circle,:circle, :circle, :circle, :circle],
            markersize = [3,5,5,5,3], 
            linewidth = 4, 
            aspect_ratio=:equal, 
            legend=false)
    end
    function PlotAnimation(x_all,y_all,title,sw)
        #x = ((x.-1).%2).+1
        n =  length(x_all[1,:])
        @gif for i ∈ 1:n
            plot([x_all[1,i]],[y_all[1,i]], markercolors=:black, markershapes=:rect, markersizes=:12)
            PlotState!(x_all[:,i],y_all[:,i],title)
            
            tot = length(sw)
            num = length(findall(sw.>i))
            if (tot - num) %2 == 0 
                plot!([x_all[5,i]-.15, x_all[5,i]+.15], [0,0], lc=:red, ls=:solid, lw=.7)
            else
                plot!([x_all[5,i]-.05, x_all[5,i]+.05], [0,0], lc=:green, ls=:solid, lw=.7)     
            end

        end every 1
    end
### END Animating ###

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
                           

##### Import Data #####
(t_all,u_all,p,filename) = ImportData()

# Generate some data
x,y = X_Y(u_all,p)
sw = Findsw(u_all)

# Animate Data
title = "Single Leg With State Transitions"
PlotAnimation(x,y,filename,sw)




