################## Sim Goal XLSX Plotter 1.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
using Makie
using GLMakie
using CairoMakie

using XLSX
using DataFrames
using LaTeXStrings

#### Plotting Parameters ####
ylab = L"\bar{V}"
Sizetitle = 5
Sizelabel = 30
Sizetick = 20
Sizeleagend = 20



### Plotting Goal Offset Data ###
# Read the Excel file
xf = XLSX.readxlsx("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\Goal_Offset_output.xlsx")

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
    xlabel = L"\alpha",
    ylabel = ylab,
    xticklabelsize = Sizetick,
    yticklabelsize = Sizetick,
    xticklabelfont = "Arial",
    yticklabelfont = "Arial",
    xlabelsize = Sizelabel,
    ylabelsize = Sizelabel,
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
    offset, vavg = data
    
    if length(vavg) == 3
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
#cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\SimPlots\\")

# Save the plot as an EPS file with the specified filename
#Makie.save("Goal_Offset_vs_AvgVelocity.eps", sp)





### Plotting Goal Px Data ###
# Read the Excel file
xf = XLSX.readxlsx("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\Goal_Px_output.xlsx")

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
    xlabel = L"P_x",
    ylabel = ylab,
    xticklabelsize = Sizetick,
    yticklabelsize = Sizetick,
    xticklabelfont = "Arial",
    yticklabelfont = "Arial",
    xlabelsize = Sizelabel,
    ylabelsize = Sizelabel,
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
    
    if length(vavg) == 3
        # Create scatter plot with specified markersize, color, and marker symbol
        scatter_plot = Makie.scatter!(ax, px, vavg, markersize = 15, color = :black, marker = symbols[i])
        
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
#cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\SimPlots\\")
cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\ExpPlots\\")
# Save the plot as an EPS file with the specified filename

Makie.save("Effect_on_AvgVelocitySimplified.png", sp)





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
    xlabelsize = Sizelabel,
    ylabelsize = Sizelabel,
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
#cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\SimPlots\\")
cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData\\ExpPlots\\")
# Save the plot as an EPS file with the specified filename
Makie.save("Test.eps", sp)