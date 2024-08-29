################## Sim Goal XLSX Plotter 1.0 ##################

################### Imported Librarys ###################
Path_name = joinpath(@__DIR__, "MyRobotFunctionPackage/src")
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end

include(joinpath(Path_name, "MyRobotFunctionPackage5.jl"))
import .MyRobotFunctionPackage5 as MF5
using Makie
using GLMakie
using CairoMakie

using XLSX
using DataFrames
using LaTeXStrings

#### Plotting Parameters ####
ylab = latexstring("Average Velocity [\$m/s\$], \$\\bar{V}\$")
Sizetitle = 5
Sizelabely = 30
Sizelabelx = 30
Sizetick = 20
Sizeleagend = 20



### Plotting Goal Offset Data ###
# Read the Excel file
parent_path = joinpath(@__DIR__, "..")
desired_folder_path = joinpath(parent_path, "Simulation Datasets")
xf = XLSX.readxlsx(joinpath(desired_folder_path, "Goal_Offset_output.xlsx"))

# Get the sheet names
sheerNames = XLSX.sheetnames(xf)

# Read the data from the first sheet
sh = xf[sheerNames[1]]
Offset = sh["A"][2:end]
Vxhipd = sh["B"][2:end]
Vavg = sh["C"][2:end]

# Separate the data into different data sets based on Vxhipd values
Vxhipd_values = unique(Vxhipd)
data_sets = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()  # Dictionary to store data sets

for vx in Vxhipd_values
    indices = findall(Vxhipd .== vx)  # Find indices where Vxhipd matches current value
    offset = Offset[indices]
    vavg = Vavg[indices]
    data_sets[vx] = (offset, vavg)  # Store data set in dictionary
end

# Create a plot
sp = Figure(resolution = (900, 700))
ax = Axis(sp[1, 1],
    title = " ", #latexstring("\$\\alpha\$ and \$P_x\$ Effect on Average Velocity"),
    titlecolor = :black,
    titlesize = Sizetitle,
    titlefont = "Times",
    xlabel = latexstring("Empirical Parameter [-], \$\\alpha\$"),
    ylabel = "",
    xticklabelsize = Sizetick,
    yticklabelsize = Sizetick,
    xticklabelfont = "Arial",
    yticklabelfont = "Arial",
    xlabelsize = Sizelabelx,
    ylabelsize = Sizelabely,
    xlabelfont = "Times",
    ylabelfont = "Times"
    )
    Label(sp[1:2, 0], ylab, 
    rotation = pi/2,
    justification = :center,
    fontsize = Sizelabely,
    font = "Times")
hidedecorations!(ax,label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

# Create empty arrays to store scatter plots and legend names
scatter_plots = []
legend_names = []

# Define symbols for different markers
symbols = [:cross, :square, :diamond, :utriangle, :dtriangle]

# Iterate over data sets
for (i, (vx, data)) in enumerate(data_sets)
    offset, vavg = data
    
    if length(vavg) >= 1
        # Create scatter plot with specified markersize, color, and marker symbol
        scatter_plot = Makie.scatter!(ax, offset, vavg, markersize = 15, color = :black, marker = symbols[i])
        
        # Add scatter plot to the array
        push!(scatter_plots, scatter_plot)
        
        # Create legend name with interpolated value of vx
        push!(legend_names, latexstring("\$\\dot{x}^{*}_{h}\$ = $(vx) m/s"))
    end
end

# Add legend to the plot with specified position
axislegend(ax, scatter_plots, legend_names, position = :rb, labelsize = Sizeleagend, labelfont = "Times")
sp

# Change current directory to the specified path
parent_path = joinpath(@__DIR__, "..\\..")
plot_folder_path = joinpath(parent_path, "Experimental Testing\\ExpPlots")
cd(plot_folder_path)






### Plotting Goal Px Data ###
# Read the Excel file
xf = XLSX.readxlsx(joinpath(desired_folder_path, "Goal_Px_output.xlsx"))

# Get the sheet names
sheerNames = XLSX.sheetnames(xf)

# Read the data from the first sheet
sh = xf[sheerNames[1]]
Px = sh["A"][2:end]
Vxhipd = sh["B"][2:end]
Vavg = sh["C"][2:end]

# Separate the data into different data sets based on Vxhipd values
Vxhipd_values = unique(Vxhipd)
data_sets = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()  # Dictionary to store data sets

for vx in Vxhipd_values
    indices = findall(Vxhipd .== vx)  # Find indices where Vxhipd matches current value
    px = Px[indices]
    vavg = Vavg[indices]
    data_sets[vx] = (px, vavg)  # Store data set in dictionary
end

# Create a plot
#sp = Figure(resolution = (1000, 500))
ax = Axis(sp[2, 1],
    title = " ", #latexstring("\$P_x\$ Effect on Average Velocity"),
    titlecolor = :black,
    titlesize = 5,
    titlefont = "Times",
    xlabel = xlabel = latexstring("Horizontal Percentage Parameter [-], \$P_x\$"),
    ylabel = "",
    xticklabelsize = Sizetick,
    yticklabelsize = Sizetick,
    xticklabelfont = "Arial",
    yticklabelfont = "Arial",
    xlabelsize = Sizelabelx,
    ylabelsize = Sizelabely,
    xlabelfont = "Times",
    ylabelfont = "Times"
    )
hidedecorations!(ax,label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

# Create empty arrays to store scatter plots and legend names
scatter_plots = []
legend_names = []

# Define symbols for different markers
symbols = [:cross, :square, :diamond, :utriangle, :dtriangle]

# Iterate over data sets
for (i, (vx, data)) in enumerate(data_sets)
    px, vavg = data
    
    if length(vavg) >= 1
        # Create scatter plot with specified markersize, color, and marker symbol
        scatter_plot = Makie.scatter!(ax, px, vavg, markersize = 15, color = :black, marker = symbols[i])
        
        # Add scatter plot to the array
        push!(scatter_plots, scatter_plot)
        
        # Create legend name with interpolated value of vx
        push!(legend_names, latexstring("\$\\dot{x}^{*}_{h}\$ = $(vx) m/s"))
    end
end

# Add legend to the plot with specified position
axislegend(ax, scatter_plots, legend_names, position = :lt, labelsize = Sizeleagend, labelfont = "Times")


sp

# Change current directory to the specified path
cd(plot_folder_path)

# Save the plot as an EPS file with the specified filename

Makie.save("Effect_on_AvgVelocity.eps", sp)




#=
sp = Figure(resolution = (1000, 500))
ax = Axis(sp[1, 1],
    title = " ", #latexstring("\$P_x\$ Effect on Average Velocity"),
    titlecolor = :black,
    titlesize = 5,
    titlefont = "Times",
    xlabel = L"P_x",
    ylabel = ylab,
    xticklabelsize = Sizetick,
    yticklabelsize = Sizetick,
    xticklabelfont = "Arial",
    yticklabelfont = "Arial",
    xlabelsize = Sizelabelx,
    ylabelsize = Sizelabely,
    xlabelfont = "Times",
    ylabelfont = "Times"
    )
hidedecorations!(ax,label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

# Create empty arrays to store line plots and legend names
line_plots = []
legend_names = []

# Define line for different plots
line_styles = [:dash, :solid, :dashdot]
line_color = [:red, :black, :black]

# Iterate over data sets
j = 0
for (i, (vx, data)) in enumerate(data_sets)
    px, vavg = data
    

    # Create line plot with specified linewidth, color, and linestyle
    line_plot = Makie.lines!(ax, px, vavg, linewidth = 2, color = line_color[i-j], linestyle = line_styles[i-j])
    
    # Add line plot to the array
    push!(line_plots, line_plot)
    
    # Create legend name with interpolated value of vx
    push!(legend_names, latexstring("\$\\dot{x}^{*}_{h}\$ = $(vx) m/s"))
end

legend_names = ["Desired ", "Measured Stance", "Measured Flight"]

# Add legend to the plot with specified position

axislegend(ax, line_plots, legend_names, position = :rb, labelsize = Sizeleagend, labelfont = "Times")
sp

# Change current directory to the specified path
cd(plot_folder_path)

# Save the plot as an EPS file with the specified filename
Makie.save("Test.eps", sp)
=#