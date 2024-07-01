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
(t_all,u_all,p,filename) = MF5.ImportData("z_DataSet99.xlsx")
# Generate some data
x1,y1 = MF5.X_Y(u_all,p)
sw1 = MF5.Findsw(u_all)

##### Import Data #####
(t_all,u_all,p,filename) = MF5.ImportData("z_DataSet999.xlsx")
# Generate some data
x2,y2 = MF5.X_Y(u_all,p)
sw2 = MF5.Findsw(u_all)

# Animate Data
title = ""
gif1 = MF5.PlotAnimation3(x1,y1,sw1,x2,y2,sw2,title)

gif(gif1, "Bounding_Simulation2.gif", fps=50)