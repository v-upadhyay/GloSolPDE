using LinearAlgebra
using StaticArrays
using Serialization
# using CairoMakie
using GLMakie

include("system.jl")

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


function exp_sin(t, a, b, terminal_time)
    integral = -b * exp(a * (terminal_time - t)) * cos(b * t) - a * exp(a * (terminal_time - t)) * sin(b * t) + b * exp(a * terminal_time)
    return integral / (a^2 + b^2)
end

function exp_cos(t, a, b, terminal_time)
    integral = -a * exp(a * (terminal_time - t)) * cos(b * t) + b * exp(a * (terminal_time - t)) * sin(b * t) + a * exp(a * terminal_time)
    return integral / (a^2 + b^2)
end


function exp_basis_integral(t, i, a, num_params, terminal_time)
    n = 10
    if i==1
        return (exp(a * (terminal_time))  - exp(a * (terminal_time - t)))/a 
    elseif i <= (num_params+1)/2
        return exp_sin(t, a, 2*pi*(i-1)/(terminal_time * n), terminal_time)
    else
        return exp_cos(t, a, 2*pi*(i-(num_params+1)/2)/(terminal_time * n), terminal_time)
    end
end

function exp_basis_span_integral(t, a, num_params, terminal_time, coeff)
    integral = 0.0
    for i in 1:num_params
        integral = integral + coeff[i] * exp_basis_integral(t, i, a, num_params, terminal_time)
    end
    return integral
end

function control(t, control_params, num_params, terminal_time)
    control_commands = 0.0
    for i in 1:num_params
        control_commands = control_commands + basis(t, i, num_params, terminal_time) * control_params[i]
    end
    return control_commands
end

function control_integral(t, control_matrix, state_eigen_values, state_eigen_vectors, control_params, num_params)
    
    d_x = length(state_eigen_values)

    integral = [exp_basis_span_integral(t, state_eigen_values[i], num_params, terminal_time, control_params) for i in 1:d_x]

    exp_integral = Diagonal(integral)
    exp_integral = state_eigen_vectors * exp_integral * inv(state_eigen_vectors)
    exp_integral = exp_integral * control_matrix
    return exp_integral
end

function uncert_integral(t, control_matrix, state_eigen_values, state_eigen_vectors, uncert_param, num_params)

    d_x = length(state_eigen_values)

    integral = [exp_basis_span_integral(t, state_eigen_values[i], num_params, terminal_time, uncert_param) for i in 1:d_x]

    exp_integral = Diagonal(integral)
    exp_integral = state_eigen_vectors * exp_integral * inv(state_eigen_vectors)
    exp_integral = exp_integral * control_matrix
    return exp_integral
end

function dynamics_constraints(initial_state, state_matrix, control_matrix, control_params, uncert, time, num_params)

    initial_distribution = initial_state

    A = state_matrix
    B = control_matrix

    psuedo_state = initial_distribution
    state = (exp(A * time)) * psuedo_state
    state  = state + control_integral(time, B, state_eigen_values, state_eigen_vectors, control_params, num_params)
    state = state + uncert_integral(time, B, state_eigen_values, state_eigen_vectors, uncert, num_params)

end

N = 50                                                                                                         # Semi-discretization grid points
state_matrix = state_params(N)                                              # State matrix
state_matrix_decomposition = eigen(state_matrix)
state_eigen_values = state_matrix_decomposition.values
state_eigen_vectors = state_matrix_decomposition.vectors
state_matrix_inv = inv(state_matrix)                                       # Inverse of the state matrix
control_matrix = control_params(N)                                          # Control matrix
num_params = 20
num_samples = num_params + 1

initial_modes = 3
terminal_mode = 1
terminal_time = 2.0                                                                                             # Terminal time
initial_state = [10 * ((x/N)^2) * (1 - x/N)^3 * sin(initial_modes * pi * x/N) for x in 1:N]  
terminal_state = [2 * ((x/N)) * (1 - x/N)^3 * sin(terminal_mode * pi * x/N) for x in 1:N]

control_1 = deserialize("Recent_Control_Trajectory.dat")

plot_samples = 100

trajectory = zeros(plot_samples, N+1)


best_sol = deserialize("disturbance_maximizer.dat")
opt_time_samples = best_sol[end-num_samples+1:end]  # Extract the last num_samples elements as time samples
opt_uncert_array = best_sol[1:end-num_samples]   # Extract the rest as uncertainty samples

# Reshape the uncertainty samples into the uncertainty tensor
opt_uncert = reshape(opt_uncert_array, num_samples, num_params)  # Reshape into uncertainty tensor

# uncert_1 = zeros(num_params)  # Nominal case without uncertainty
uncert_1 = opt_uncert[1, :]  # Extract the first row as the uncertainty sample

# Load parameters and compute control
control_params_array = control_1 + uncert_1  # Use the first row of the uncertainty tensor as control parameters
t = LinRange(0, terminal_time, plot_samples)
control_array = [control(ti, control_params_array, num_params, terminal_time) for ti in t]


for i in 1:plot_samples
    trajectory[i, :] = vcat(dynamics_constraints(initial_state, state_matrix, control_matrix, control_1, uncert_1, i * terminal_time / plot_samples, num_params), [control_array[i]])
end

# Plotting the trajectory

t_axis = LinRange(0, terminal_time, plot_samples)
x_axis = LinRange(0, 1, N+1)

terminal_controlled = zeros(N)
for i in 1:N
    terminal_controlled[i] = trajectory[plot_samples, i]
end

# Plot using CairoMakie
fig = Figure(size = (900, 640), backgroundcolor = :white)
label_size = 28
tick_size = 16
ax = Axis3(fig[1, 1],
    xlabel = "Space", ylabel = "Time", zlabel = "Amplitude",
    xlabelsize = label_size, ylabelsize = label_size, zlabelsize = label_size,
    xticklabelsize = tick_size, yticklabelsize = tick_size, zticklabelsize = tick_size,
    xlabelcolor = :black, ylabelcolor = :black, zlabelcolor = :black,
    xgridvisible = true, ygridvisible = true, zgridvisible = true,
    xgridcolor = (:gray, 0.3), ygridcolor = (:gray, 0.3), zgridcolor = (:gray, 0.3),
    xgridwidth = 1, ygridwidth = 1, zgridwidth = 1,
    backgroundcolor = :white,
    title = "Evolution of State over Time",
    titlesize = 20
)
surface!(ax, x_axis, t_axis, trajectory', colormap = :inferno, shading = true)
fig

save("state_trajectory.png", fig)
fig

# fig = Figure()
# ax = Axis(fig[1, 1])

# lines!(ax, terminal_controlled; label = "Terminal Controlled")
# lines!(ax, terminal_state; label = "Terminal State")

# axislegend(ax)
# fig
