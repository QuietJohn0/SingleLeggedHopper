################## Double Pendulum XLSX Animator 2.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
using Plots
using ProgressMeter

##### Import Data #####
(t_all,u_all,p,filename) = MF5.ImportData()

# Generate some data
x,y = MF5.X_Y(u_all,p)
sw = MF5.Findsw(u_all)

# Animate Data
title = "Bounding Simulation"
gif1 = MF5.PlotAnimation(x,y,title,sw)

gif(gif1, "Bounding_Simulation.gif", fps=50)