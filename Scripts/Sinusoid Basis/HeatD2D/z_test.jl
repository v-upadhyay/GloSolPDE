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

# Load nominal control and disturbance data
nominal_control_params = deserialize("Recent_Control_Trajectory.dat")
opt_uncert = deserialize("disturbance_maximizer.dat")

# Isolate the specific disturbance vector
disturbance_params = opt_uncert[1, :]

# Create the disturbed control parameters by adding the disturbance
disturbed_control_params = nominal_control_params + disturbance_params

# --- Compute Both Control Trajectories for Comparison ---

# Generate time steps
t = LinRange(0, terminal_time, 1000)

# Calculate the control values for both scenarios
nominal_control = [control(ti, nominal_control_params, num_params, time_period) for ti in t]
disturbed_control = [control(ti, disturbed_control_params, num_params, time_period) for ti in t]

# Load nominal control and disturbance data
nominal_control_params_1 = deserialize("Recent_Control_Trajectory_1.dat")
opt_uncert_1 = deserialize("disturbance_maximizer_1.dat")

# Isolate the specific disturbance vector
disturbance_params_1 = opt_uncert_1[1, :]

# Create the disturbed control parameters by adding the disturbance
disturbed_control_params_1 = nominal_control_params_1 + disturbance_params_1

# Calculate the control values for both scenarios
nominal_control_1 = [control(ti, nominal_control_params_1, num_params, time_period) for ti in t]
disturbed_control_1 = [control(ti, disturbed_control_params_1, num_params, time_period) for ti in t]

# --- Enhanced Plotting ---

# Create a figure with a slightly larger size for better clarity
fig = Figure(resolution = (900, 600))
label_size = 30
tick_size = 18
# Add an axis with more descriptive labels and title
ax = Axis(fig[1, 1],
    xlabel = "Time (s)",
    ylabel = "Boundary Input",
    xlabelsize = label_size, ylabelsize = label_size, 
    xticklabelsize = tick_size, yticklabelsize = tick_size, 
    titlesize = 30,

    # NEW: Enable and style the X and Y grids
    xgridvisible = true,
    ygridvisible = true,
    xgridstyle = :dash,
    ygridstyle = :dash,
    xgridcolor = RGBA(0, 0, 0, 0.15), # Light, semi-transparent black
    ygridcolor = RGBA(0, 0, 0, 0.15),

    # Keep the modern look by hiding top/right borders
    topspinevisible = false,
    rightspinevisible = false,

    # NEW: Increase tick density for a finer grid
    xticks = WilkinsonTicks(10),
    yticks = WilkinsonTicks(10)
)

# # Plot the nominal trajectory (dashed gray line)
# lines!(ax, t, nominal_control,
#     color = :,
#     linewidth = 2.5,
#     linestyle = :dash,
#     label = "Nominal Control"
# )

# Plot the disturbed trajectory (solid blue line)
# lines!(ax, t, disturbed_control,
#     color = :dodgerblue,
#     linewidth = 3.5,
#     label = "Disturbed Control"
# )

lines!(ax, t, disturbed_control_1,
    color = :gray,
    linewidth = 3.5,
    label = "Disturbed Control_1"
)

# # Add a shaded band to emphasize the difference
# band!(ax, t, nominal_control, disturbed_control,
#     color = RGBA(0.2, 0.6, 1.0, 0.25),
#     label = "Disturbance Effect"
# )

# Add a vertical line to mark the terminal time
vlines!(ax, [terminal_time],
    color = :black,
    linestyle = :dot,
    linewidth = 2,
    label = "Terminal Time"
)

# NEW: Customize the legend with a larger font size
axislegend(ax,
    position = :rb,
    framevisible = false,
    labelsize = 20  # Increased from default
)

# Optional: Tighten axis limits to the data range for a snug fit
xlims!(ax, 0, terminal_time)

# Save the enhanced figure with higher resolution for sharpness
save("HeatContD2D.png", fig, px_per_unit = 2)

# Display the figure
# fig