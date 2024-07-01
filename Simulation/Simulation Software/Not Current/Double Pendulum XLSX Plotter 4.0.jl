################## Double Pendulum XLSX Plotter 3.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
using Plots

##### State Transition Plots #####
    function PlotTransitionLines(t_all,u, title,sw)
        plot(t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], title=title, legend=false)
        umin = findmin(u[sw[1]+1:sw[end]])[1]
        umax = findmax(u[sw[1]+1:sw[end]])[1]
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                plot!([t_all[j], t_all[j]], [umin, umax], lc=:red, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [umin, umax], lc=:red, ls=:dash, lw=.7)
            elseif j != k
                plot!([t_all[j], t_all[j]], [umin, umax], lc=:green, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [umin, umax], lc=:green, ls=:dash, lw=.7)        
            end
        end
        plot!([t_all[sw[end]] t_all[sw[end]]], [umin umax], lc=:green, lw=1)
    end
    function PlotTransitionLines!(t_all,u, title,sw)
        plot!(t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], title=title, legend=false)
        umin = findmin(u[sw[1]+1:sw[end]])[1]
        umax = findmax(u[sw[1]+1:sw[end]])[1]
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                plot!([t_all[j], t_all[j]], [umin, umax], lc=:red, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [umin, umax], lc=:red, ls=:dash, lw=.7)
            elseif j != k
                plot!([t_all[j], t_all[j]], [umin, umax], lc=:green, ls=:dash, lw=.7)
                plot!([t_all[k-1], t_all[k-1]], [umin, umax], lc=:green, ls=:dash, lw=.7)        
            end
        end
        plot!([t_all[sw[end]] t_all[sw[end]]], [umin umax], lc=:green, lw=1)
    end
    function PlotTransitionPoints(ux,uy, title,sw)
        plot(ux[sw[1]+1:sw[end]], uy[sw[1]+1:sw[end]], title=title, legend=false)
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
        plot!([ux[sw[end]] ux[sw[end]]], [uy[sw[end]] uy[sw[end]]], markercolors = :black, shape = :circle, markersize = 2)
    end
    function PlotTransitionPoints!(ux,uy, title,sw)
        plot!(ux[sw[1]+1:sw[end]], uy[sw[1]+1:sw[end]], title=title, legend=false)
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
        plot!([ux[sw[end]] ux[sw[end]]], [uy[sw[end]] uy[sw[end]]], markercolors = :black, shape = :circle, markersize = 2)
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
        PlotTransitionLines(t_all, -p.kₛ*u_all[:,10] - p.cₛ*u[:,12], title, sw)
        Fd = p.A*sin.(p.w.*t_all+ones(size(t_all))*p.ϕ)
        PlotTransitionLines!(t_all, Fd, title, sw)
    end
    function PlotStancePy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,4], title, sw)
        plot!([t_all[1], t_all[end]], p.yhipd*[1,1])
    end
    function PlotStanceVx(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,7], title, sw)
        plot!([t_all[1], t_all[end]], p.Vxhipd*[1,1])
    end
    function PlotStanceVy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,8], title, sw)
        plot!([t_all[1], t_all[end]], p.Vyhipd*[1,1])
    end
### END Stance Controller Plots ###

##### Stance Controller Plots #####
function PlotAngularMomentum(t_all,p,u_all,sw,title)
    H = MF5.Calculate_Angular_Momentum_All(p,u_all,sw,true)
    PlotTransitionLines(t_all, H[1,:], title, sw)
    PlotTransitionLines!(t_all, H[2,:], title, sw)
end
function TrimData(sw,Start_Cycle,N_Cycles)
    ## Flight then Stance 1 cycle ##
    sw[Start_Cycle*2+1:(Start_Cycle+N_Cycles)*2+2]
end
function PlotCycleAngularMomentum(t_all,p,u_all,sw,title)
    sw = TrimData(sw,0,1)
    PlotAngularMomentum(t_all,p,u_all,sw,title)
end
### END Stance Controller Plots ###



##### Import Data #####
(t_all,u_all,p) = MF5.ImportData()

# Generate some data
x,y = MF5.X_Y(u_all,p)
sw = MF5.Findsw(u_all)
τ = MF5.findTorques(u_all,p,t_all,sw)
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
p4 = PlotTransitionLines(t_all, p.kₛ*u_all[:,10], "grf", sw)
p5 = PlotTransitionPoints(u_all[:,1], u_all[:,1] + u_all[:,2], "θ₁ v.s. θ₁+θ₃", sw)
p6 = PlotTransitionPoints(t_all, τ[1,:], "Torque, τ", sw)
p6 = PlotTransitionPoints!(t_all, τ[2,:], "Torque, τ", sw)
sp2 = plot(p1, p2, p3, p4, p5, p6, layout=(2,3))
#display(sp2)

##### Subplot 3: Stance Controls #####
p1 = PlotStanceGRF(t_all,u_all,p, "GRF",sw)
p2 = PlotStancePy(t_all,u_all,p, "Pyhid",sw)
p3 = PlotStanceVx(t_all,u_all,p, "Vxhid",sw)
p4 = PlotStanceVy(t_all,u_all,p, "Vyhid",sw)
sp3 = plot(p1, p2, p3, p4, layout=(2,2))
display(p3)

##### Subplot 4: Conservation of Angular Momentum #####
#p1 = PlotAngularMomentum(t_all,p,u_all,sw,"Angular Momentum")
#p2 = PlotCycleAngularMomentum(t_all,p,u_all,sw,"Angular Momentum 2 Cycles")

#p1 = PlotSegmentCompare(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)
