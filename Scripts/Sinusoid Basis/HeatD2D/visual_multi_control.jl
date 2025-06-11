# Import necessary libraries
using GLMakie
using LinearAlgebra
using Serialization
using Colors

# Apply a clean, modern theme for a professional look
# This theme customizes fonts, removes the top/right axis spines, and styles the grid.
set_theme!(theme_minimal())

# Parameters (unchanged)
terminal_time = 10.0
num_params = 20
time_period = 10 * terminal_time

# Basis function (unchanged)
function basis(t, i, num_params, time_period)
    if i == 1
        return 1
    elseif i <= (num_params + 1) / 2
        return sin(2 * pi * (i - 1) * t / time_period)
    else
        return cos(2 * pi * (i - (num_params + 1) / 2) * t / time_period)
    end
end

# Control computation (unchanged)
function control(t, control_params, num_params, time_period)
    control_commands = 0.0
    for i in 1:num_params
        control_commands += basis(t, i, num_params, time_period) * control_params[i]
    end
    return control_commands
end

# --- Data Loading and Processing ---



num_params = 20 
num_samples = num_params + 1

using ColorSchemes
num_lines = 21  # Number of lines to plot
# Number of samples
num_colors = num_lines


# --- Enhanced Plotting ---

N_sample = [20, 50, 100, 150, 200]

for N in N_sample


# Load nominal control and disturbance data
nominal_control_params = deserialize("Recent_Control_Trajectory_$N.dat")
opt_uncert = deserialize("disturbance_maximizer_$N.dat")

# Choose a color scheme
cmap = get(ColorSchemes.coolwarm, range(0, 1, length=num_colors))

# Create a figure with enhanced resolution and styling
fig = Figure(resolution = (900, 600))
label_size = 30
tick_size = 18

ax = Axis(fig[1, 1],
    xlabel = "Time (s)",
    ylabel = "Boundary Input",
    xlabelsize = label_size, ylabelsize = label_size, 
    xticklabelsize = tick_size, yticklabelsize = tick_size, 
    titlesize = 30,
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash,
    xgridcolor = RGBA(0, 0, 0, 0.15),
    ygridcolor = RGBA(0, 0, 0, 0.15),
    topspinevisible = false,
    rightspinevisible = false,
    xticks = WilkinsonTicks(10),
    yticks = WilkinsonTicks(10)
)

# Plot disturbed inputs with gradient colors
for i in 1:num_lines
    disturbance_params = opt_uncert[i, :]
    disturbed_control_params = nominal_control_params + disturbance_params
    t = LinRange(0, terminal_time, 1000)
    disturbed_control = [control(ti, disturbed_control_params, num_params, time_period) for ti in t]

    lines!(ax, t, disturbed_control,
        color = cmap[i],
        linewidth = 3.5,
        label = "Input_$i"
    )
end

# Mark terminal time
vlines!(ax, [terminal_time],
    color = :black,
    linestyle = :dot,
    linewidth = 2,
    label = "Terminal Time"
)

xlims!(ax, 0, terminal_time)
save("HeatContD2D_$N.png", fig, px_per_unit = 2)

fig

end