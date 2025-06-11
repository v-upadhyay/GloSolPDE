using Ipopt
using LinearAlgebra
using StaticArrays

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


function control_constraints(model, control_params, uncert, num_params, control_bounds, time_samples, time_period)
    
    lower_bound = control_bounds[1]
    upper_bound = control_bounds[2]

    for i in 1:size(time_samples)[1]

        for j in 1:size(time_samples)[1]
        @constraint(model, lower_bound <= control(time_samples[i], control_params + uncert[j, :], num_params, time_period) <= upper_bound)
        end

    end

    @constraint(model, control(0, control_params, num_params, time_period) + control(0, uncert[1, :], num_params, time_period) == 0.0)

end

function slack_constraint(model, terminal_state, lambda, gamma, control_params, slack, states_sampled, uncert, num_params, terminal_time)
    
    for i in 1:num_samples
    
        states_sampled_error = states_sampled[:, i] - terminal_state
        @constraint(model, lambda * sum(control_params.^2) + 0.5 * 100 * sum(states_sampled_error.^2) - gamma * sum(uncert[i, :].^2) <=  slack)

    end
end

function dynamics_constraints(model, initial_state, state_matrix, control_matrix, states_sampled, control_params, uncert, num_samples, num_params, time_period, terminal_time)

    initial_distribution = initial_state

    A = state_matrix
    B = control_matrix
    
    for i in 1:num_samples

    psuedo_state = initial_distribution
    state = (exp(A * terminal_time)) * psuedo_state
    state  = state + control_integral(terminal_time, B, state_eigen_values, state_eigen_vectors, control_params, num_params, time_period)
    state = state + uncert_integral(terminal_time, B, state_eigen_values, state_eigen_vectors, uncert[i,:], num_params, time_period)

    @constraint(model, states_sampled[:, i] .== state)

    end
end