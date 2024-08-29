################## Double Pendulum XLSX Animator 2.0 ##################

################### Imported Librarys ###################
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
x1,y1 = MF5.X_Y(u_all,p)
sw1 = MF5.Findsw(u_all)

##### Import Data #####
(t_all,u_all,p,filename) = MF5.ImportData("CurrentDatasetTuned2.xlsx")
# Generate some data
x2,y2 = MF5.X_Y(u_all,p)
sw2 = MF5.Findsw(u_all)

# Animate Data
title = ""
gif1 = MF5.PlotAnimation3(x1,y1,sw1,x2,y2,sw2,title)

parent_path = joinpath(@__DIR__, "..")
desired_folder_path = joinpath(parent_path, "Animations")
cd(desired_folder_path)

gif(gif1, "Bounding_Simulation2.gif", fps=50)