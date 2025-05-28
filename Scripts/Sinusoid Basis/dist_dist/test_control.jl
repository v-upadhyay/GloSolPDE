# using CairoMakie
using GLMakie
using LinearAlgebra
using Serialization

# Parameters
terminal_time = 2
num_params = 20
num_samples = num_params + 1
# Basis function
function basis(t, i, num_params, terminal_time)
    n = 10
    if i==1
        return 1
    elseif i <= (num_params+1)/2
        return sin(2*pi*(i-1)*t/(terminal_time * n))
    else
        return cos(2*pi*(i-(num_params+1)/2)*t/(terminal_time * n))
    end
end

# Control computation
function control(t, control_params, num_params)
    u = 0.0
    for i in 1:num_params
        u += basis(t, i, num_params, terminal_time) * control_params[i]
    end
    return u
end

opt_uncert = deserialize("disturbance_maximizer.dat")


uncert_1 = opt_uncert[1, :]  # Extract the first row as the uncertainty sample

# Load parameters and compute control
control_params_array = deserialize("Recent_Control_Trajectory.dat") + uncert_1  # Use the first row of the uncertainty tensor as control parameters
# control_params_array = uncert_1
# control_params_array = opt_uncert[2, :]  # Use the first row of the uncertainty tensor as control parameters
t = LinRange(0, terminal_time, 1000)
control_array = [control(ti, control_params_array, num_params) for ti in t]

# Plot using CairoMakie
fig = Figure(resolution = (800, 500), backgroundcolor = :white)
ax = Axis(fig[1, 1],
    xlabel = "Time", ylabel = "Control Value",
    xlabelsize = 24, ylabelsize = 24,
    xticklabelsize = 18, yticklabelsize = 18,
    title = "Control Trajectory", titlesize = 26,
    xgridvisible = true, ygridvisible = true,
    xgridcolor = (:gray, 0.3), ygridcolor = (:gray, 0.3)
)

lines!(ax, t, control_array, color = :dodgerblue, linewidth = 3, label = "Control")
axislegend(ax, position = :rt, textsize = 18, labelsize = 18)

save("control_trajectory.png", fig)
fig
