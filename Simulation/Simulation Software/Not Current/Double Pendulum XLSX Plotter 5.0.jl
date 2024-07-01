################## Double Pendulum XLSX Plotter 5.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
using Plots
using LaTeXStrings
using MathJaxRenderer


using ColorSchemes 
using Plots # Get the Viridis color scheme 
#viridis_colors = ColorSchemes.viridis(10) # Plot using the Viridis color scheme 
#plot(rand(10), rand(10), color=viridis_colors, seriestype=:scatter)

##### State Transition Plots #####
    function PlotSetUp(t_all, u, sw, title, xlab, ylab)
        umin = findmin(u[sw[1]+1:sw[end]])[1]
        umax = findmax(u[sw[1]+1:sw[end]])[1]
        Gap = (umax - umin)*.05

        if title == " "
            titlesize = 5
        else
            titlesize = 20
        end

        plot(xlims = (t_all[sw[1]+1], t_all[sw[end]]), 
             ylims = (umin-Gap, umax+Gap), 
             grid=false,
             framestyle = :box,
             legend=false,
             title = title,
             titlefont=font(titlesize,"Times"),
             xlabel = xlab,
             ylabel = ylab,
             tickfont = font(10,"Times"),
             guidefont=font(12,"Arial"),
             extra_plot_kwargs = KW(:include_mathjax => "cdn"),
             dpi = 300
             )
    end
    function PlotSetUp(t_all, u1, u2, sw, title, xlab, ylab)
        umin = findmin([u1[sw[1]+1:sw[end]]; u2])[1]
        umax = findmax([u1[sw[1]+1:sw[end]]; u2])[1]       

        PlotSetUp([t_all[sw[1]+1] t_all[sw[end]]], [umin umax], [0 2], title, xlab, ylab)
    end
    function TransitionLines!(t_all,u,sw)
        umin = findmin(u[sw[1]+1:sw[end]])[1]
        umax = findmax(u[sw[1]+1:sw[end]])[1]
        Gap = (umax - umin)*.05

        linestyle = :dash
        linecolor = :grey
        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i == 1
                plot!([t_all[j+1], t_all[j+1]], [umin-Gap, umax+Gap], lc=linecolor, ls=linestyle, lw=.7)
                plot!([t_all[k], t_all[k]], [umin-Gap, umax+Gap], lc=linecolor, ls=linestyle, lw=.7)
            elseif i%2 == 0
                plot!([t_all[j], t_all[j]], [umin-Gap, umax+Gap], lc=linecolor, ls=linestyle, lw=.7)
                plot!([t_all[k], t_all[k]], [umin-Gap, umax+Gap], lc=linecolor, ls=:dash, lw=.7)
            elseif j != k
                plot!([t_all[j], t_all[j]], [umin-Gap, umax+Gap], lc=linecolor, ls=linestyle, lw=.7)
                plot!([t_all[k], t_all[k]], [umin-Gap, umax+Gap], lc=linecolor, ls=linestyle, lw=.7)        
            end
        end
        plot!([t_all[sw[end]] t_all[sw[end]]], [umin umax], lc=linecolor, lw=1)
    end
    function PlotTransitionLines(t_all, u, sw, title, xlab, ylab)
        PlotSetUp(t_all, u, sw, title, xlab, ylab)
        plot!(t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], lc = :black, lw = 1)
        TransitionLines!(t_all,u,sw)
    end
    function PlotTransitionLines!(t_all, u, sw)
        plot!(t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], lc = :black, lw = 1)
        TransitionLines!(t_all,u,sw)
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
    function PlotTransitionColor(t_all, u, sw, title, xlab, ylab)
        PlotSetUp(t_all, u, sw, title, xlab, ylab)

        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i == 1
                plot!(t_all[j+1:k], u[j+1:k], lc=:black, ls =:dashdot ,lw = 1)
            elseif i%2 == 0
                plot!(t_all[j:k], u[j:k], lc=:black, ls =:solid, lw = 2)
            elseif j != k
                plot!(t_all[j:k], u[j:k], lc=:black, ls =:dashdot, lw = 1)   
            end
        end
        plot!()
    end
    function PlotTransitionColor!(t_all, u, sw)
        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i == 1
                plot!(t_all[j+1:k], u[j+1:k], lc=:black, ls =:dashdot ,lw = 1)
            elseif i%2 == 0
                plot!(t_all[j:k], u[j:k], lc=:black, ls =:solid, lw = 2)
            elseif j != k
                plot!(t_all[j:k], u[j:k], lc=:black, ls =:dashdot, lw = 1)   
            end
        end
        plot!()
    end
    function PlotTransitionLinesColor(t_all, u, sw, title, xlab, ylab)
        PlotTransitionColor(t_all, u, sw, title, xlab, ylab)
        TransitionLines!(t_all,u,sw)
    end
### END State Transition Plots ###

##### Stance Controller Plots #####
    function PlotStanceGRF(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, -p.kₛ*u_all[:,10] - p.cₛ*u_all[:,12].*(u_all[:,12].<0), sw, "",L"t, \text{[s]}", L"GRF, \text{[N]}")
        Fd = p.A*sin.(p.w.*t_all+ones(size(t_all))*p.ϕ)
        PlotTransitionLines!(t_all, Fd, sw)
    end
    function PlotStancePy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,4], title, sw)
        plot!([t_all[1], t_all[end]], p.yhipd*[1,1])
    end
    function PlotVx!(t_all, u, sw, p)
        plot!([t_all[sw[1]], t_all[sw[end]]], p.Vxhipd*[1,1])
    end
    function PlotStanceVy(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, u_all[:,8], title, sw)
        plot!([t_all[1], t_all[end]], p.Vyhipd*[1,1])
    end
### END Stance Controller Plots ###

##### Stance Controller Plots #####
    function PlotAngularMomentum(t_all, u_all, sw, p, title, xlab, ylab)
        H = MF5.Calculate_Angular_Momentum_All(p,u_all,sw,true)
        p1 = PlotTransitionLinesColor(t_all, H[1,:], sw, title, "", L"P_x, \text{[N-s]}")
        p2 = PlotTransitionLinesColor(t_all, H[2,:], sw, " ", "", L"P_y, \text{[N-s]}")
        p3 = PlotTransitionLinesColor(t_all, H[3,:], sw, "", "", L"H_A, \text{[N-m-s]}")
        p4 = PlotTransitionLinesColor(t_all, H[4,:], sw, " ", xlab, L"H_B, \text{[N-m-s]}")
        plot(p1, p2, p3, p4, layout=(4,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
    end
    function TrimData(sw,Start_Cycle,N_Cycles)
        ## Flight then Stance is 1 cycle ##
        sw[round(Int32,Start_Cycle*2+1):round(Int32,(Start_Cycle+N_Cycles)*2+2)]
    end
    function PlotCycleAngularMomentum(t_all, u_all, sw, p, title, xlab, ylab)
        sw = TrimData(sw,0,1)
        PlotAngularMomentum(t_all, u_all, sw, p, title, xlab, ylab)
    end
### END Stance Controller Plots ###


##### Import Data #####
(t_all,u_all,p) = MF5.ImportData()

# Generate some data
x,y = MF5.X_Y(u_all,p)
sw = MF5.Findsw(u_all)
τ = MF5.findTorques(u_all,p,t_all,sw)
GRF = -p.kₛ*u_all[:,10] - p.cₛ*u_all[:,12].*(u_all[:,12].<0)
plotlyjs()

# Plotting Parameters
swplot = TrimData(sw,0,8.5)

##### Subplot 1: x, y #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,3], swplot, "", "", L"x_0, \text{[m]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,4], swplot, " ", L"t, \text{[s]}", L"y_0, \text{[m]}")
sp1 = plot(p1, p2, layout=(2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
#savefig("test.eps")

##### Subplot 2: ẋ, ẏ #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,7], swplot, "" ,""         ,L"\dot{x}_0, \text{[m/s]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,9], swplot, " ",L"t, \text{[s]}",L"\dot{y}_0, \text{[m/s]}")
sp2 = plot(p1, p2, layout=(2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))

##### Subplot 3: θ₁, θ₃, θ̇₁, θ̇₃ #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,1], swplot, "" ,"",L"\theta_1, \text{[rad]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,2], swplot, " ","",L"\theta_2, \text{[rad]}")
p3 = PlotTransitionLinesColor(t_all, u_all[:,5], swplot, " ","",L"\dot{\theta}_1, \text{[rad/s]}")
p4 = PlotTransitionLinesColor(t_all, u_all[:,6], swplot, " ",L"t, \text{[s]}",L"\dot{\theta}_2, \text{[rad/s]}")
sp3 = plot(p1, p2, p3, p4, layout=(4,1),size = (700,500),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
#savefig("test.html")

##### Subplot 4: GRF #####
sp4 = PlotTransitionLinesColor(t_all, GRF, swplot, "",L"t, \text{[s]}", L"GRF, \text{[N]}")

##### Subplot 5: τ₁, τ₂ #####
p1 = PlotTransitionLinesColor(t_all, τ[1,:], swplot, "" ,"", L"\tau_1, \text{[N/m]}")
p2 = PlotTransitionLinesColor(t_all, τ[2,:], swplot, " ",L"t, \text{[s]}", L"\tau_2, \text{[N/m]}")
sp5 = plot(p1, p2, layout = (2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))

##### Subplot 6/7: CON MOM #####
sp6 = PlotCycleAngularMomentum(t_all, u_all, sw, p, "", L"t, \text{[s]}", "")
sp7 = PlotAngularMomentum(t_all, u_all, swplot, p, "", L"t, \text{[s]}", "")
savefig("AngularMomentumPlot.html")

##### Display #####
for i in [sp1, sp2,sp3,sp4,sp5,sp6,sp7]; display(i); end

##### Subplot 8: Flight #####
swflight = TrimData(sw,24,0)
(t, θ₁d, θ₃d, θ̇₁d, θ̇₃d) = MF5.Calculate_Desired_Flight(p,u_all[swflight[1],:],t_all[swflight[1],:])

p1 = PlotSetUp(t_all, u_all[:,1], θ₁d, swflight, "" ,"",L"\theta_1, \text{[rad]}")
p1 = plot!(t, θ₁d)
p1 = PlotTransitionColor!(t_all, u_all[:,1], swflight)

#less y labels
p2 = PlotSetUp(t_all, u_all[:,2],  θ₃d, swflight, " ","",L"\theta_2, \text{[rad]}")
p2 = plot!(t, θ₃d)
p2 = yticks!([-.5, -.7, -.9, -1.1])
p2 = ylims!((-1.1, -.5))
p2 = PlotTransitionColor!(t_all, u_all[:,2], swflight)

p3 = PlotSetUp(t_all, u_all[:,5], θ̇₁d, swflight, " ","",L"\dot{\theta}_1, \text{[rad/s]}")
p3 = plot!(t, θ̇₁d)
p3 = PlotTransitionColor!(t_all, u_all[:,5], swflight)

p4 = PlotSetUp(t_all, u_all[:,6], θ̇₃d, swflight, " ",L"t, \text{[s]}",L"\dot{\theta}_2, \text{[rad/s]}")
p4 = plot!(t, θ̇₃d)
p4 = PlotTransitionColor!(t_all, u_all[:,6], swflight)

sp8 = plot(p1, p2, p3, p4, layout=(4,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
plot!(size = (700,500))
#savefig("test.html")

##### Subplot 9/10: Stance #####
swstance = TrimData(sw,24.5,0)
(t, θ₁d, θ₃d, θ̇₁d, θ̇₃d, GRFd) = MF5.Calculate_Desired_Stance(p,u_all[swstance[1],:],t_all[swstance[1],:])

p1 = PlotSetUp(t_all, u_all[:,1], θ₁d, swstance, "" ,"",L"\theta_1, \text{[rad]}")
p1 = plot!(t, θ₁d)
p1 = PlotTransitionColor!(t_all, u_all[:,1], swstance)

p2 = PlotSetUp(t_all, u_all[:,2], θ₃d, swstance, " ","",L"\theta_2, \text{[rad]}")
p2 = plot!(t, θ₃d)
p2 = PlotTransitionColor!(t_all, u_all[:,2], swstance)

p3 = PlotSetUp(t_all, u_all[:,5], θ̇₁d, swstance, " ","",L"\dot{\theta}_1, \text{[rad/s]}")
p3 = plot!(t, θ̇₁d)
p3 = PlotTransitionColor!(t_all, u_all[:,5], swstance)

p4 = PlotSetUp(t_all, u_all[:,6], θ̇₃d, swstance, " ",L"t, \text{[s]}",L"\dot{\theta}_2, \text{[rad/s]}")
p4 = plot!(t, θ̇₃d)
p4 = PlotTransitionColor!(t_all, u_all[:,6], swstance)

sp9 = plot(p1, p2, p3, p4, layout=(4,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
plot!(size = (700,500))

sp10 = PlotSetUp(t_all, GRF, GRFd, swstance, " ",L"t, \text{[s]}",L"GRF_d, \text{[N]}")
sp10 = plot!(t, GRFd)
sp10 = PlotTransitionColor!(t_all, GRF, swstance)




# Steady State Plotting Parameters
swplot = TrimData(sw,30,2)

##### Subplot 11: x, y #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,3], swplot, "", "", L"x_0, \text{[m]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,4], swplot, " ", L"t, \text{[s]}", L"y_0, \text{[m]}")
sp11 = plot(p1, p2, layout=(2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))

##### Subplot 12: ẋ, ẏ #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,7], swplot, "" ,""         ,L"\dot{x}_0, \text{[m/s]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,8], swplot, " ",L"t, \text{[s]}",L"\dot{y}_0, \text{[m/s]}")
sp12 = plot(p1, p2, layout=(2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))


p1 = PlotSetUp(t_all, u_all[:,7], p.Vxhipd, swplot, "" ,""         ,L"\dot{x}_0, \text{[m/s]}")
p1 = PlotVx!(t_all, u_all[:,7], swplot,p)
p1 = PlotTransitionColor!(t_all, u_all[:,7], swplot)


##### Subplot 13: θ₁, θ₃, θ̇₁, θ̇₃ #####
p1 = PlotTransitionLinesColor(t_all, u_all[:,1], swplot, "" ,"",L"\theta_1, \text{[rad]}")
p2 = PlotTransitionLinesColor(t_all, u_all[:,2], swplot, " ","",L"\theta_2, \text{[rad]}")
p3 = PlotTransitionLinesColor(t_all, u_all[:,5], swplot, " ","",L"\dot{\theta}_1, \text{[rad/s]}")
p4 = PlotTransitionLinesColor(t_all, u_all[:,6], swplot, " ",L"t, \text{[s]}",L"\dot{\theta}_2, \text{[rad/s]}")
sp13 = plot(p1, p2, p3, p4, layout=(4,1),size = (700,500),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
#savefig("test.html")

##### Subplot 14: GRF #####
sp14 = PlotTransitionLinesColor(t_all, GRF, swplot, "",L"t, \text{[s]}", L"GRF, \text{[N]}")

##### Subplot 15: τ₁, τ₂ #####
p1 = PlotTransitionLinesColor(t_all, τ[1,:], swplot, "" ,"", L"\tau_1, \text{[N/m]}")
p2 = PlotTransitionLinesColor(t_all, τ[2,:], swplot, " ",L"t, \text{[s]}", L"\tau_2, \text{[N/m]}")
sp15 = plot(p1, p2, layout = (2,1),size = (700,500),extra_plot_kwargs = KW(:include_mathjax => "cdn"))


##### Subplot 16: Trajectory Planning #####
    swflight = TrimData(sw,24,0)
    (t, θ₁d, θ₃d, θ̇₁d, θ̇₃d) = MF5.Calculate_Desired_Flight(p,u_all[swflight[1],:],t_all[swflight[1],:])
    t = t .- t[1]

    p1 = PlotSetUp(t_all, u_all[:,1], θ₁d, swflight, "" ,"",L"\theta_{1d}, \text{[rad]}")
    p1 = plot!(t, θ₁d, lc=:black, ls =:dashdot, lw = 1)
    p1 = scatter!([t[1], t[end]], [θ₁d[1], θ₁d[end]], markercolors = :red, shape = :circle, markersize = 2)


    p3 = PlotSetUp(t_all, u_all[:,5], θ̇₁d, swflight, " ","",L"\dot{\theta}_{1d}, \text{[rad/s]}")
    p3 = plot!(t, θ̇₁d, lc=:black, ls =:dashdot, lw = 1)
    p3 = scatter!([t[1], t[end]], [θ̇₁d[1], θ̇₁d[end]], markercolors = :red, shape = :circle, markersize = 2)


    sp8 = plot(p1, p3, layout=(2,1),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
    sp8 = xlims!((t[1],t[end]))
    plot!(size = (700,500))
    #savefig("test.html")
##### Subplot 16: Plotted #####



















##### Subplot 1: State Space #####
p1 = PlotTransitionLines(t_all, u_all[:,1], "θ₁", sw)
p2 = PlotTransitionLines(t_all, u_all[:,2], "θ₃", sw)
p3 = PlotTransitionLines(t_all, u_all[:,3], "x", swplot)
p4 = PlotTransitionLines(t_all, u_all[:,4], "y", swplot)
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
p6 = PlotTransitionPoints(t_all, τ[2,:], "Torque, τ", sw)
sp2 = plot(p1, p2, p3, p4, p5, p6, layout=(2,3))
#display(sp2)

##### Subplot 3: Stance Controls #####
p1 = PlotStanceGRF(t_all,u_all,p, "GRF",swstance)
p2 = PlotStancePy(t_all,u_all,p, "Pyhid",sw)
p3 = PlotStanceVx(t_all,u_all,p, "Vxhid",sw)
p4 = PlotStanceVy(t_all,u_all,p, "Vyhid",sw)
sp3 = plot(p1, p2, p3, p4, layout=(2,2))
display(p3)

##### Subplot 4: Conservation of Angular Momentum #####
#p1 = PlotAngularMomentum(t_all,p,u_all,sw,"Angular Momentum")
#p2 = PlotCycleAngularMomentum(t_all,p,u_all,sw,"Angular Momentum 2 Cycles")

#p1 = PlotSegmentCompare(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)
