################## Double Pendulum XLSX Plotter 6.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
#using Plots
using Makie
using GLMakie
using CairoMakie


#using ColorSchemes 
#using Plots # Get the Viridis color scheme 
#viridis_colors = ColorSchemes.viridis(10) # Plot using the Viridis color scheme 
#plot(rand(10), rand(10), color=viridis_colors, seriestype=:scatter)

##### State Transition Plots #####
    function PlotSetUp(fig, t_all, u, sw, title, xlab, ylab)
        if sw[1] == 0
            sw[1] = 1
        end
        
        umin = findmin(u[sw[1]:sw[end]])[1]
        umax = findmax(u[sw[1]:sw[end]])[1]
        Gap = (umax - umin) * 0.05

        if title == " "
            titlesize = 5
        else
            titlesize = 20
        end

        ax = Axis(fig,
            title = title,
            titlecolor = :black,
            titlesize = titlesize*2,
            titlefont = "Times",
            xlabel = xlab,
            ylabel = ylab,
            limits = (t_all[sw[1]], t_all[sw[end]],umin - Gap, umax + Gap),
            xticklabelsize = 10*2,
            yticklabelsize = 10*2,
            xticklabelfont = "Arial",
            yticklabelfont = "Arial",
            xlabelsize = 16*2,
            ylabelsize = 16*2,
            xlabelfont = "Times",
            ylabelfont = "Times"
        )

        hidedecorations!(ax,label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

        ax
    end
    function PlotSetUp(fig,t_all, u1, u2, sw, title, xlab, ylab)
        umin = findmin([u1[sw[1]+1:sw[end]]; u2])[1]
        umax = findmax([u1[sw[1]+1:sw[end]]; u2])[1]       

        ax = PlotSetUp(fig,[t_all[sw[1]+1] t_all[sw[end]]], [umin umax], [0 2], title, xlab, ylab)
    end
    function TransitionLines!(ax,t_all,u,sw)
        if sw[1] == 0
            sw[1] = 1
        end
        umin = findmin(u[sw[1]:sw[end]])[1]
        umax = findmax(u[sw[1]:sw[end]])[1]
        Gap = (umax - umin)*.05

        lstyle = :dash
        linecolor = :grey
        
        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i == 1
                lines!(ax, [t_all[j], t_all[j]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2)
                lines!(ax, [t_all[k], t_all[k]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2)
            elseif i%2 == 0
                #lines!(ax, [t_all[j], t_all[j]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
                lines!(ax, [t_all[k], t_all[k]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2)
            elseif j != k
                #lines!(ax, [t_all[j], t_all[j]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
                lines!(ax, [t_all[k], t_all[k]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2) 
            end
        end
        #lines!(ax, [t_all[sw[end]], t_all[sw[end]]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
        
    end
    function PlotTransitionLines(fig,t_all, u, sw, title, xlab, ylab)
        ax = PlotSetUp(fig,t_all, u, sw, title, xlab, ylab)
        lines!(ax, t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], color = :black, linewidth = 2)
        TransitionLines!(t_all,u,sw)
        ax
    end
    function PlotTransitionLines!(ax,t_all, u, sw)
        lines!(ax, t_all[sw[1]+1:sw[end]], u[sw[1]+1:sw[end]], color = :black, linewidth = 2)
        TransitionLines!(ax,t_all,u,sw)
    end
    function PlotTransitionColor(fig,t_all, u, sw, title, xlab, ylab)
        ax = PlotSetUp(fig,t_all, u, sw, title, xlab, ylab)
        TransitionColor!(ax,t_all, u, sw)
        ax
    end
    #=function TransitionColor!(ax,t_all, u, sw)
        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i == 1
                lines!(ax, t_all[j:k], u[j:k], color = :black, linestyle = :dashdot, linewidth = 2)
                #plot!(t_all[j+1:k], u[j+1:k], lc=:black, ls =:dashdot ,lw = 1)
            elseif i%2 == 0
                lines!(ax, t_all[j:k], u[j:k], color = :black, linestyle = :solid, linewidth = 3)
                #plot!(t_all[j:k], u[j:k], lc=:black, ls =:solid, lw = 2)
            elseif j != k
                lines!(ax, t_all[j:k], u[j:k], color = :black, linestyle = :dashdot, linewidth = 2)
                #plot!(t_all[j:k], u[j:k], lc=:black, ls =:dashdot, lw = 1)   
            end
        end
    end=#
    function TransitionColor!(ax,t_all, u, sw, linecolor = :black, lstyle = :dashdot, lwidth = 2)  
        if sw[1] == 0
            sw[1] = 1
        end
        if lstyle == :solid
            style = :dashdot
            width = 2
        elseif lstyle == :dashdot
            style = :solid
            width = 3
        else
            style = lstyle
            width = lwidth
        end

        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i%2 == 0
                lines!(ax, t_all[j:k], u[j:k], color = linecolor, linestyle = style, linewidth = width)                #plot!(t_all[j:k], u[j:k], lc=:black, ls =:dashdot, lw = 1)
            elseif j != k
                lines!(ax, t_all[j:k], u[j:k], color = linecolor, linestyle = lstyle, linewidth = lwidth)
                #plot!(t_all[j:k], u[j:k], lc=:black, ls =:solid, lw = 2)   
            end
        end
    end
    function PlotTransitionLinesColor(fig,t_all, u, sw, title, xlab, ylab)
        ax = PlotTransitionColor(fig,t_all, u, sw, title, xlab, ylab)
        TransitionLines!(ax,t_all,u,sw)
        ax
    end
### END State Transition Plots ###

##### Old State Transition Plots #####
    function PlotSetUpOld(t_all, u, sw, title, xlab, ylab)
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
            dpi = 900,
            )
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
### END Old State Transition Plots ###

##### Stance Controller Plots #####
    function PlotStanceGRF(t_all,u,p, title,sw)
        PlotTransitionLines(t_all, -p.kₛ*u_all[:,10] - p.cₛ*u_all[:,12].*(u_all[:,12].<0), sw, "",L"t \text{, [s]}", L"GRF \text{, [N]}")
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
    function PlotAngularMomentum(t_all, u_all, sw, swplot, p, title, xlab, ylab)
        H = MF5.Calculate_Angular_Momentum_All(p,u_all,sw,true)
        sp = Figure(resolution = (900, 600))
        
        #p1 = PlotTransitionLinesColor(sp[1,1],t_all, H[1,:], swplot, title, "", L"P_x \text{, [N-s]}")
        #p1 = yticks!([3.4, 3.3, 3.2, 3.1])
        #p1 = ylims!((3.1,3.4))
        #p2 = PlotTransitionLinesColor(sp[2,1],t_all, H[2,:], swplot, " ", "", L"P_y \text{, [N-s]}")
        #p2 = yticks!([1.5, .5, -.5, -1.5])
        p3 = PlotTransitionLinesColor(sp[1,1],t_all, H[3,:], swplot, "", "", L"H_A \text{, [Nms]}")
        p4 = PlotTransitionLinesColor(sp[2,1],t_all, H[4,:], swplot, " ", xlab, L"H_B \text{, [Nms]}")
        #p4 = yticks!([-.75, -.8, -.85])

        #plot(p1, p2, p3, p4, layout=(4,1),size = (700,500),extra_plot_kwargs = KW(:include_mathjax => "cdn"))
        sp
    end
    function TrimData(sw,Start_Cycle,N_Cycles)
        ## Flight then Stance is 1 cycle ##
        sw[round(Int32,Start_Cycle*2+1):round(Int32,(Start_Cycle+N_Cycles)*2+2)]
    end
    function PlotCycleAngularMomentum(t_all, u_all, sw, p, title, xlab, ylab)
        swplot = TrimData(sw,20,1)
        PlotAngularMomentum(t_all, u_all, sw, swplot, p, title, xlab, ylab)
    end
### END Stance Controller Plots ###







##### Import Data #####
(t_all,u_all,p) = MF5.ImportData()

(value, index) = findmax(u_all[:,4])
#index = 1
t_all = t_all[index:end] .- t_all[index]
u_all = u_all[index:end,:]

# Generate some data
x,y = MF5.X_Y(u_all,p)
sw = MF5.Findsw(u_all)
τ = MF5.findTorques(u_all,p,t_all,sw)
GRF = (-p.kₛ*u_all[:,10] - p.cₛ*u_all[:,12].*(u_all[:,12].<0)).*(round.(u_all[:,10]*1E15)/1E15.<=0)

# Plotting Parameters
swplot = TrimData(sw,0,9.5)

cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\SimPlots\\")
##### Subplot 1: x, y #####
sp1 = Figure(resolution = (900, 700))
PlotTransitionLinesColor(sp1[1,1],t_all, u_all[:,3], swplot, "", "", L"x_0 \text{, [m]}")
PlotTransitionLinesColor(sp1[2,1],t_all, u_all[:,4], swplot, " ", L"t \text{, [s]}", L"y_0 \text{, [m]}")
sp1
Makie.save("HipPosition_Trans.eps", sp1)

##### Subplot 2: ẋ, ẏ #####
sp2 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp2[1,1],t_all, u_all[:,7], swplot, "" ,""         ,L"\dot{x}_0 \text{, [m/s]}")
PlotTransitionLinesColor(sp2[2,1],t_all, u_all[:,8], swplot, " ",L"t \text{, [s]}",L"\dot{y}_0 \text{, [m/s]}")
sp2
Makie.save("HipVelocity_Trans.eps", sp2)

##### Subplot 3: θ₁, θ₃, θ̇₁, θ̇₃ #####
sp3 = Figure(resolution = (900, 800))
PlotTransitionLinesColor(sp3[1,1], t_all, u_all[:,1], swplot, "" ,"",L"\theta_1 \text{, [rad]}")
PlotTransitionLinesColor(sp3[2,1], t_all, u_all[:,2], swplot, " ","",L"\theta_2 \text{, [rad]}")
PlotTransitionLinesColor(sp3[3,1], t_all, u_all[:,5], swplot, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
PlotTransitionLinesColor(sp3[4,1], t_all, u_all[:,6], swplot, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
sp3
Makie.save("AngularPosition_Trans.eps", sp3)

##### Subplot 4: GRF #####
sp4 = Figure(resolution = (900, 400))
PlotTransitionLinesColor(sp4[1,1], t_all, GRF, swplot, "",L"t \text{, [s]}", L"R_{y}, \text{ [N]}")
sp4
Makie.save("GRF_Trans.eps", sp4)

##### Subplot 5: τ₁, τ₂ #####
sp5 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp5[1,1], t_all, τ[1,:], swplot, "" ,"", L"\tau_1 \text{, [Nm]}")
PlotTransitionLinesColor(sp5[2,1], t_all, τ[2,:], swplot, " ",L"t \text{, [s]}", L"\tau_2 \text{, [Nm]}")
sp5
Makie.save("Torque_Trans.eps", sp5)

###### Subplot 6/7: CON MOM ######
sp6 = PlotCycleAngularMomentum(t_all, u_all, sw, p, "", L"t \text{, [s]}", "")
#sp7 = PlotAngularMomentum(t_all, u_all, swplot, p, "", L"t \text{, [s]}", "")
Makie.save("AngularMomentumPlot.eps", sp6)
##### End Subplot 6/7: CON MOM #####


##### Subplot 8: Flight #####
swflight = TrimData(sw,24,0)
(t, θ₁d, θ₃d, θ̇₁d, θ̇₃d) = MF5.Calculate_Desired_Flight(p,u_all[swflight[1],:],t_all[swflight[1],:],swflight[end]-swflight[1])

sp8 = Figure(resolution = (900, 800))

ax = PlotSetUp(sp8[1,1],t_all, u_all[:,1], θ₁d, swflight, "" ,"",L"\theta_1 \text{, [rad]}")
lines!(ax, t, θ₁d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,1], swflight)

ax = PlotSetUp(sp8[2,1], t_all, u_all[:,2],  θ₃d, swflight, " ","",L"\theta_2 \text{, [rad]}")
lines!(ax, t, θ₃d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,2], swflight)

ax = PlotSetUp(sp8[3,1], t_all, u_all[:,5], θ̇₁d, swflight, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
lines!(ax, t, θ̇₁d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,5], swflight)

ax = PlotSetUp(sp8[4,1], t_all, u_all[:,6], θ̇₃d, swflight, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
lines!(ax, t, θ̇₃d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,6], swflight)
sp8

Makie.save("Flight_SS.eps", sp8)

##### Subplot 9/10: Stance #####
swstance = TrimData(sw,24.5,0)
(t, θ₁d, θ₃d, θ̇₁d, θ̇₃d, GRFd) = MF5.Calculate_Desired_Stance(p,u_all[swstance[1],:],t_all[swstance[1],:])

sp9 = Figure(resolution = (900, 800))

ax = PlotSetUp(sp9[1,1],t_all, u_all[:,1], θ₁d, swstance, "" ,"",L"\theta_1 \text{, [rad]}")
lines!(ax, t, θ₁d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,1], swstance, :black, :solid, 2)

ax = PlotSetUp(sp9[2,1], t_all, u_all[:,2],  θ₃d, swstance, " ","",L"\theta_2 \text{, [rad]}")
lines!(ax, t, θ₃d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,2], swstance, :black, :solid, 2)

ax = PlotSetUp(sp9[3,1], t_all, u_all[:,5], θ̇₁d, swstance, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
lines!(ax, t, θ̇₁d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,5], swstance, :black, :solid, 2)

ax = PlotSetUp(sp9[4,1], t_all, u_all[:,6], θ̇₃d, swstance, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
lines!(ax, t, θ̇₃d, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, u_all[:,6], swstance, :black, :solid, 2)
sp9
Makie.save("Stance_SS.eps", sp9)

sp10 = Figure(resolution = (900, 400))
swstance[1] = swstance[1] - 1
ax = PlotSetUp(sp10[1,1], t_all, GRF, GRFd, swstance, " ",L"t \text{, [s]}",L"R_{y_d} \text{, [N]}")
lines!(ax, t, GRFd, color = :red, linewidth = 2, linestyle = :dash)
TransitionColor!(ax, t_all, GRF, swstance, :black, :solid, 2)
sp10
Makie.save("Stance_SSGRF.eps", sp10)



# Steady State Plotting Parameters
swplotss = TrimData(sw,24,1)

##### Subplot 11: x, y #####
sp11 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp11[1,1], t_all, u_all[:,3], swplotss, "", "", L"x_0 \text{, [m]}")
PlotTransitionLinesColor(sp11[2,1], t_all, u_all[:,4], swplotss, " ", L"t \text{, [s]}", L"y_0 \text{, [m]}")
sp11
Makie.save("HipPosition_SS.eps", sp11)

##### Subplot 12: ẋ, ẏ #####
sp12 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp12[1,1], t_all, u_all[:,7], swplotss, "" ,""         ,L"\dot{x}_0 \text{, [m/s]}")
PlotTransitionLinesColor(sp12[2,1], t_all, u_all[:,8], swplotss, " ",L"t \text{, [s]}",L"\dot{y}_0 \text{, [m/s]}")
sp12
Makie.save("HipVelocity_SS.eps", sp12)

#p1 = PlotSetUp(t_all, u_all[:,7], p.Vxhipd, swplotss, "" ,""         ,L"\dot{x}_0 \text{, [m/s]}")
#p1 = PlotVx!(t_all, u_all[:,7], swplotss,p)
#p1 = TransitionColor!(t_all, u_all[:,7], swplotss)


##### Subplot 13: θ₁, θ₃, θ̇₁, θ̇₃ #####
sp13 = Figure(resolution = (900, 800))
PlotTransitionLinesColor(sp13[1,1], t_all, u_all[:,1], swplotss, "" ,"",L"\theta_1 \text{, [rad]}")
PlotTransitionLinesColor(sp13[2,1], t_all, u_all[:,2], swplotss, " ","",L"\theta_2 \text{, [rad]}")
PlotTransitionLinesColor(sp13[3,1], t_all, u_all[:,5], swplotss, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
PlotTransitionLinesColor(sp13[4,1], t_all, u_all[:,6], swplotss, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
sp13
Makie.save("AngularPosition_SS.eps", sp13)

##### Subplot 14: GRF #####
sp14 = Figure(resolution = (900, 400))
PlotTransitionLinesColor(sp14[1,1], t_all, GRF, swplotss, "",L"t \text{, [s]}", L"R_{y} \text{, [N]}")
sp14
Makie.save("GRF_SS.eps", sp14)

##### Subplot 15: τ₁, τ₂ #####
sp15 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp15[1,1], t_all, τ[1,:], swplotss, "" ,"", L"\tau_1 \text{, [Nm]}")
PlotTransitionLinesColor(sp15[2,1], t_all, τ[2,:], swplotss, " ",L"t \text{, [s]}", L"\tau_2 \text{, [Nm]}")
sp15
Makie.save("Torque_SS.eps", sp15)

##### Subplot 16: Trajectory Planning #####
    swflight = TrimData(sw,24,0)
    (t, θ₁d, θ₃d, θ̇₁d, θ̇₃d) = MF5.Calculate_Desired_Flight(p,u_all[swflight[1],:],t_all[swflight[1],:],swflight[end]-swflight[1])
    t = t .- t[1]

    sp16 = Figure(resolution = (900, 600))

    ax = PlotSetUp(sp16[1,1], t, θ₁d, swflight .- swflight[1], "" ,"",L"\theta^{*} \text{, [rad]}")
    lines!(ax, t, θ₁d, color = :black, linestyle =:dashdot, linewidth = 3)
    Makie.scatter!(ax, [t[1], t[end]], [θ₁d[1], θ₁d[end]], color = :red, markersize = 20)

    ax = PlotSetUp(sp16[2,1], t, θ̇₁d, swflight .- swflight[1], "" ,"",L"\dot{\theta}^{*} \text{, [rad]}")   
    lines!(ax, t, θ̇₁d, color = :black, linestyle =:dashdot, linewidth = 3)
    Makie.scatter!(ax, [t[1], t[end]], [θ̇₁d[1], θ̇₁d[end]], color = :red, markersize = 20)
    sp16
    Makie.save("TragectoryPlanningExample.eps", sp16)
##### Subplot 16: Plotted #####

















#=

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

#p1 = PlotSegmentCompare(u_all[:,1],u_all[:,5], "θ₁ v.s. θ̇₁",sw)=#
