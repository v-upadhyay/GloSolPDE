using LinearAlgebra
using StaticArrays
using Serialization
# using CairoMakie
using GLMakie

include("system.jl")

function basis(t, i, num_params, time_period)
    if i==1
        return 1
    elseif i <= (num_params+1)/2
        return sin(2*pi*(i-1)*t/time_period)
    else
        return cos(2*pi*(i-(num_params+1)/2)*t/time_period)
    end
end


function exp_sin(t, a, b)
    integral = -b * cos(b * t) - a * sin(b * t) + b * exp(a * t)
    return integral / (a^2 + b^2)
end

function exp_cos(t, a, b)
    integral = - a * cos(b * t) + b * sin(b * t) + a * exp(a * t)
    return integral / (a^2 + b^2)
end


function exp_basis_integral(t, i, a, num_params, time_period)
    # t: current integration limit
    # a: eigenvalue
    # terminal_time_overall_sim: This argument from the original signature is now unused
    #                              within this function's logic due to the correction.
    #                              It was the source of the issue.
    # time_period_basis: period for the sin/cos basis functions

    if i == 1
        # We need Integral_0^t exp(a*(t-s))*1 ds
        # This evaluates to (exp(a*t) - 1)/a for a != 0, and t for a = 0.
        if isapprox(a, 0.0, atol=1e-9) # Handle a=0 case to avoid division by zero
            return t
        else
            return (exp(a * t) - 1.0) / a
        end
    elseif i <= (div(num_params + 1, 2)) # Integer division for clarity
        # We need Integral_0^t exp(a*(t-s))*sin(freq*s) ds
        # The original exp_sin(t_limit, alpha, beta_freq, T_in_exp_factor)
        # calculates Integral_0^{t_limit} exp(alpha*(T_in_exp_factor - s))*sin(beta_freq*s)ds.
        # To get the correct form, T_in_exp_factor must be t_limit.
        # Here, t_limit is the first argument 't'.
        freq = 2 * pi * (i - 1) / (time_period)
        return exp_sin(t, a, freq) # Pass 't' as the 4th argument to exp_sin
    else
        # Similarly for exp_cos
        freq = 2 * pi * (i - div(num_params + 1, 2)) / (time_period)
        return exp_cos(t, a, freq) # Pass 't' as the 4th argument to exp_cos
    end
end

function exp_basis_span_integral(t, a, num_params, coeff, time_period)
    integral = 0.0
    for i in 1:num_params
        integral = integral + coeff[i] * exp_basis_integral(t, i, a, num_params, time_period)
    end
    return integral
end

function control(t, control_params, num_params, time_period)
    control_commands = 0.0
    for i in 1:num_params
        control_commands = control_commands + basis(t, i, num_params, time_period) * control_params[i]
    end
    return control_commands
end

function control_integral(t, control_matrix, state_eigen_values, state_eigen_vectors, control_params, num_params, time_period)
    
    d_x = length(state_eigen_values)

    integral = [exp_basis_span_integral(t, state_eigen_values[i], num_params, control_params, time_period) for i in 1:d_x]

    exp_integral = Diagonal(integral)
    exp_integral = state_eigen_vectors * exp_integral * inv(state_eigen_vectors)
    exp_integral = exp_integral * control_matrix
    return exp_integral
end

function uncert_integral(t, control_matrix, state_eigen_values, state_eigen_vectors, uncert_param, num_params, time_period)

    d_x = length(state_eigen_values)

    integral = [exp_basis_span_integral(t, state_eigen_values[i], num_params, uncert_param, time_period) for i in 1:d_x]

    exp_integral = Diagonal(integral)
    exp_integral = state_eigen_vectors * exp_integral * inv(state_eigen_vectors)
    exp_integral = exp_integral * control_matrix
    return exp_integral
end

function dynamics_constraints(initial_state, state_matrix, control_matrix, control_params, uncert, time, num_params, time_period)

    initial_distribution = initial_state

    A = state_matrix
    B = control_matrix

    psuedo_state = initial_distribution
    state = (exp(A * time)) * psuedo_state
    state  = state + control_integral(time, B, state_eigen_values, state_eigen_vectors, control_params, num_params, time_period)
    state = state + uncert_integral(time, B, state_eigen_values, state_eigen_vectors, uncert, num_params, time_period)

end


N_og = 100  
N_samples = 15  # Number of samples to test
error_N = zeros(N_samples)  # Initialize error_N to store terminal state errors for each N
N_plot = [i * N for i in 1:N_samples]  # N values for plotting

num_params = 20
num_samples = num_params + 1

initial_modes = 3
terminal_mode = 1
terminal_time = 0.3                                                                                             # Terminal time
time_period = 10 * terminal_time                                                                                 # Time period for the control


control_1 = deserialize("Recent_Control_Trajectory.dat")
    
    
opt_uncert = deserialize("disturbance_maximizer.dat")

# uncert_1 = zeros(num_params)  # Nominal case without uncertainty
uncert_1 = opt_uncert[1, :]  # Extract the first row as the uncertainty sample

for i in 1:N_samples

    N = N_plot[i]  # Update N for each sample

    state_matrix = state_params(N)                                              # State matrix
    control_matrix = control_params(N)                                          # Control matrix

    terminal_state = [0.5 * (x/N) for x in 1:N]                                                                            # Initial distribution for the PDE
    initial_state = initial_state = [10 * ((x/N)^2) * (1 - x/N)^3 * sin(initial_modes * pi * x/N) for x in 1:N]                                                     # Terminal distribution for the PDE
    
    error_N[i] = norm(dynamics_constraints(initial_state, state_matrix, control_matrix, control_1, uncert_1, terminal_time, num_params, time_period) - terminal_state)  # Calculate the error between controlled terminal state and desired terminal state

end

                                                                                                       # Semi-discretization grid points



# fig = Figure()
# ax = Axis(fig[1, 1])


# lines!(ax, error_N; label = "Actual Terminal State")

# axislegend(ax)
# fig
