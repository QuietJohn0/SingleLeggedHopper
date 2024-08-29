################## Double Pendulum XLSX Animator 2.0 ##################

################### Imported Librarys ###################
## Change the current working directory to the desired folder
Path_name = joinpath(@__DIR__, "MyRobotFunctionPackage/src")
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end

include(joinpath(Path_name, "MyRobotFunctionPackage5.jl"))
import .MyRobotFunctionPackage5 as MF5
using Plots
using ProgressMeter

##### Import Data #####
parent_path = joinpath(@__DIR__, "..")
desired_folder_path = joinpath(parent_path, "Simulation Datasets")
cd(desired_folder_path)

(t_all,u_all,p,filename) = MF5.ImportData("CurrentDatasetTuned1.xlsx")

# Generate some data
x,y = MF5.X_Y(u_all,p)
sw = MF5.Findsw(u_all)

# Animate Data
title = "Bounding Simulation"
gif1 = MF5.PlotAnimation(x,y,title,sw)

parent_path = joinpath(@__DIR__, "..")
desired_folder_path = joinpath(parent_path, "Animations")
cd(desired_folder_path)

gif(gif1, "Bounding_Simulation.gif", fps=50)