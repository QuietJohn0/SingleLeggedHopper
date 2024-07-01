################## Double Pendulum XLSX EXP Plotter 2.0 ##################

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
        umin = findmin([u1[sw[1]+1:sw[end]]; u2[sw[1]+1:sw[end]]])[1]
        umax = findmax([u1[sw[1]+1:sw[end]]; u2[sw[1]+1:sw[end]]])[1]       

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
            if i%2 == 0
                #lines!(ax, [t_all[j], t_all[j]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
                lines!(ax, [t_all[k], t_all[k]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2)
            elseif j != k
                #lines!(ax, [t_all[j], t_all[j]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
                lines!(ax, [t_all[k], t_all[k]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=1.2) 
            end
        end
        #lines!(ax, [t_all[sw[end]], t_all[sw[end]]], [umin-Gap, umax+Gap], color=linecolor, linestyle=lstyle, linewidth=.7)
        
    end
    function TransitionLines!(ax,t_all,u1,u2,sw)
        if sw[1] == 0
            sw[1] = 1
        end
        umin = findmin([u1[sw[1]+1:sw[end]]; u2[sw[1]+1:sw[end]]])[1]
        umax = findmax([u1[sw[1]+1:sw[end]]; u2[sw[1]+1:sw[end]]])[1] 
        Gap = (umax - umin)*.05

        lstyle = :dash
        linecolor = :grey

        for i in 1:length(sw)-1
            j = sw[i]
            k = sw[i+1]
            if i%2 == 0
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
    function PlotTransitionColor(fig,t_all, u, sw, title, xlab, ylab, linecolor = :black, lstyle = :dashdot, lwidth = 2)
        ax = PlotSetUp(fig,t_all, u, sw, title, xlab, ylab)
        TransitionColor!(ax,t_all, u, sw, linecolor, lstyle, lwidth)
        ax
    end
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
        p3 = PlotTransitionLinesColor(sp[1,1],t_all, H[3,:], swplot, "", "", L"H_A \text{, [N-m-s]}")
        p4 = PlotTransitionLinesColor(sp[2,1],t_all, H[4,:], swplot, " ", xlab, L"H_B \text{, [N-m-s]}")
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
# 8-11 40 consecutive hops
(t_all,u_all,p,τ_ff,I,θd,Id,sw,state) = MF5.ImportExpData(9)


# Generate some data
x,y = MF5.X_Y(u_all,p)
(τ, GRF, GRFc) = MF5.MotorCurrentInterpret(p,I,u_all,sw)

GRFd = zeros(size(GRF))
for i in 1:length(GRF[:,1])
    GRFd[i,:] = MF5.InvFeedForwardTorque(p,u_all[i,:],τ_ff[i,:])
end
GRFd[GRFd .< 0] .= 0

avg_u = u -> (u[1:end-4] + u[2:end-3] + u[3:end-2] + u[4:end-1] + u[5:end])./5
avg_t = t -> t[1:end-4]


# Plotting Parameters
swplot = TrimData(sw,0.5,9.5)

# Steady State Plotting Parameters
swplotss = TrimData(sw,15.5,1)


cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\ExpPlots\\")
##### Subplot 1: x, y #####
sp1 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp1[1,1],avg_t(t_all), avg_u(u_all[:,3]), swplot, "", "", L"x_0 \text{, [m]}")
PlotTransitionLinesColor(sp1[2,1],avg_t(t_all), avg_u(u_all[:,4]), swplot, " ", L"t \text{, [s]}", L"y_0 \text{, [m]}")
sp1
Makie.save("HipPosition_Trans.png", sp1)

##### Subplot 3: θ₁, θ₃, θ̇₁, θ̇₃ #####
sp3 = Figure(resolution = (900, 800))
PlotTransitionLinesColor(sp3[1,1], avg_t(t_all), avg_u(u_all[:,1]), swplot, "" ,"",L"\theta_1 \text{, [rad]}")
PlotTransitionLinesColor(sp3[2,1], avg_t(t_all), avg_u(u_all[:,2]), swplot, " ","",L"\theta_2 \text{, [rad]}")
PlotTransitionLinesColor(sp3[3,1], avg_t(t_all), avg_u(u_all[:,5]), swplot, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
PlotTransitionLinesColor(sp3[4,1], avg_t(t_all), avg_u(u_all[:,6]), swplot, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
sp3
Makie.save("AngularPosition_Trans.eps", sp3)

##### Subplot 4: GRF #####
sp4 = Figure(resolution = (900, 400))
ax = PlotTransitionLinesColor(sp4[1,1], t_all, GRF, swplot, "",L"t \text{, [s]}", L"R_{y}, \text{ [N]}")
ax.yticks = [0, 2, 4, 6]
sp4
Makie.save("GRF_Trans.png", sp4)

##### Subplot 5: I₁, I₂ #####
sp5 = Figure(resolution = (900, 600))
PlotTransitionLinesColor(sp5[1,1], avg_t(t_all), avg_u(I[:,1]), swplot, "" ,"", L"I_1 \text{, [Nm]}")
PlotTransitionLinesColor(sp5[2,1], avg_t(t_all), avg_u(I[:,2]), swplot, " ",L"t \text{, [s]}", L"I_2 \text{, [Nm]}")
sp5
Makie.save("Current_Trans.eps", sp5)





##### Subplot 11: x, y #####
    sp11 = Figure(resolution = (900, 600))
    PlotTransitionLinesColor(sp11[1,1], t_all, u_all[:,3], swplotss, "", "", L"x_0 \text{, [m]}")
    PlotTransitionLinesColor(sp11[2,1], t_all, u_all[:,4], swplotss, " ", L"t \text{, [s]}", L"y_0 \text{, [m]}")
    sp11

    Makie.save("HipPosition_SS.png", sp11)
    
##### Subplot 11: x, y #####


#p1 = PlotSetUp(t_all, u_all[:,7], p.Vxhipd, swplotss, "" ,""         ,L"\dot{x}_0 \text{, [m/s]}")
#p1 = PlotVx!(t_all, u_all[:,7], swplotss,p)
#p1 = TransitionColor!(t_all, u_all[:,7], swplotss)


##### Subplot 13: θ₁, θ₃, θ̇₁, θ̇₃ #####
sp13 = Figure(resolution = (900, 800))
PlotTransitionLinesColor(sp13[1,1], avg_t(t_all), avg_u(u_all[:,1]), swplotss, "" ,"",L"\theta_1 \text{, [rad]}")
PlotTransitionLinesColor(sp13[2,1], avg_t(t_all), avg_u(u_all[:,2]), swplotss, " ","",L"\theta_2 \text{, [rad]}")
PlotTransitionLinesColor(sp13[3,1], avg_t(t_all), avg_u(u_all[:,5]), swplotss, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
PlotTransitionLinesColor(sp13[4,1], avg_t(t_all), avg_u(u_all[:,6]), swplotss, " ",L"t \text{, [s]}",L"\dot{\theta}_2 \text{, [rad/s]}")
sp13
Makie.save("AngularPosition_SS.eps", sp13)

##### Subplot 14: GRF #####
sp14 = Figure(resolution = (900, 600))

ax = PlotSetUp(sp14[1:3,1],t_all, GRFd[:,2], GRF, swplotss, "","", L"R_{y} \text{, [N]}")
TransitionLines!(sp14[1:3,1], t_all, GRFd[:,2], swplotss)
TransitionColor!(sp14[1:3,1], t_all, GRFd[:,2], swplotss, :red, :dash, 2)
TransitionColor!(sp14[1:3,1], t_all, GRF, swplotss)

ax = PlotTransitionLinesColor(sp14[4,1], t_all, state[:,1], swplotss, " " ,L"t \text{, [s]}", L"state")
ax.yticks = [-2, 0, 2]
sp14
Makie.save("GRF_SS.eps", sp14)

##### Subplot 16: state I₁, I₂ #####
sp16 = Figure(resolution = (900, 700))
PlotTransitionLinesColor(sp16[1:2,1], t_all, I[:,1], swplotss, "" ,"", L"I_1 \text{, [A]}")
PlotTransitionLinesColor(sp16[3:4,1], t_all, I[:,2], swplotss, " ","", L"I_2 \text{, [A]}")
ax = PlotTransitionLinesColor(sp16[5,1], t_all, state[:,1], swplotss, " " ,L"t \text{, [s]}", L"state")
ax.yticks = [-2, 0, 2]
sp16
Makie.save("Current_SS.eps", sp16)

##### Subplot 17: İ, state, ẏ #####
n = 3
ẏ = u_all[1,8]
for i in u_all[1:end,8]
    ẏ = [ẏ; ẏ[end]*(n-1)/n+i/n]
end
sp17 = Figure(resolution = (900, 800))
PlotTransitionLinesColor(sp17[1:3,1], t_all, Id[:,1], swplotss, " ","", L"\dot{I}_1 \text{, [A/s]}")
PlotTransitionLinesColor(sp17[4:6,1], t_all, Id[:,2], swplotss, " ","", L"\dot{I}_2 \text{, [A/s]}")
PlotTransitionLinesColor(sp17[7:10,1], t_all, ẏ, swplotss, " ","", L"\dot{y}_0 \text{, [m/s]}")
ax = PlotTransitionLinesColor(sp17[11:12,1], t_all, state[:,1], swplotss, "" ,L"t \text{, [s]}", L"state")
ax.yticks = [-2, 0, 2]
sp17
Makie.save("StateTransition_SS.eps", sp17)

##### Subplot 8: Desired Position #####
sp8 = Figure(resolution = (900, 1100))

ax = PlotSetUp(sp8[1:2,1],t_all, u_all[:,1], θd[:,1], swplotss, "" ,"",L"\theta_1 \text{, [rad]}")
TransitionLines!(ax,t_all, u_all[:,1], θd[:,1], swplotss)    
TransitionColor!(ax, t_all, u_all[:,1], swplotss,:black,:dashdot,2)
TransitionColor!(ax, t_all, θd[:,1], swplotss,:red,:dash,2)

ax = PlotSetUp(sp8[3:4,1], t_all, u_all[:,2],  θd[:,2], swplotss, " ","",L"\theta_2 \text{, [rad]}")
TransitionLines!(ax,t_all, u_all[:,2], θd[:,2], swplotss)    
TransitionColor!(ax, t_all, u_all[:,2], swplotss,:black,:dashdot,2)
TransitionColor!(ax, t_all, θd[:,2], swplotss,:red,:dash,2)

ax = PlotSetUp(sp8[5:6,1], t_all, u_all[:,5], θd[:,3], swplotss, " ","",L"\dot{\theta}_1 \text{, [rad/s]}")
TransitionLines!(ax,t_all, u_all[:,5], θd[:,4], swplotss)
TransitionColor!(ax, t_all, u_all[:,5], swplotss,:black,:dashdot,2)
TransitionColor!(ax, t_all, θd[:,3], swplotss,:red,:dash,2)

ax = PlotSetUp(sp8[7:8,1], t_all, u_all[:,6], θd[:,4], swplotss, " ","",L"\dot{\theta}_2 \text{, [rad/s]}")
TransitionLines!(ax,t_all, u_all[:,6], θd[:,4], swplotss)
TransitionColor!(ax, t_all, u_all[:,6], swplotss,:black,:dashdot,2)
TransitionColor!(ax, t_all, θd[:,4], swplotss,:red,:dash,2)

ax = PlotTransitionLinesColor(sp8[9,1], t_all, state[:,1], swplotss, " " ,L"t \text{, [s]}", L"state")
ax.yticks = [-2, 0, 2]

ax = PlotTransitionLinesColor(sp8[9,1], t_all, state[:,1], swplotss, " " ,L"t \text{, [s]}", L"state")
ax.yticks = [-2, 0, 2]
legend!(sp8[9,1], ["State"])

sp8

Makie.save("Desired_Position_SS.eps", sp8)
## END Desired Position ##
