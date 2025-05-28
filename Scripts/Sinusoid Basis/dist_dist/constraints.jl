using Ipopt
using LinearAlgebra
using StaticArrays

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
    integral = - a * exp(a * (terminal_time - t)) * cos(b * t) + b * exp(a * (terminal_time - t)) * sin(b * t) + a * exp(a * terminal_time)
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


function control_constraints(model, control_params, num_params, control_bounds, time_samples)
    
    lower_bound = control_bounds[1]
    upper_bound = control_bounds[2]

    for i in 1:size(time_samples)[1]

        @constraint(model, lower_bound <= control(time_samples[i], control_params, num_params, terminal_time) <= upper_bound)
        # @constraint(model, control(0, control_params, num_params, terminal_time) == 0.0)

    end

end


# function slack_constraint(model, terminal_state, lambda, gamma, control_params, slack, states_sampled, uncert, num_params, terminal_time)
    
#     for i in 1:num_samples
    
#         states_sampled_error = states_sampled[:, i] - terminal_state
#         @constraint(model, lambda * sum(control_params.^2) + 0.5 * sum(states_sampled_error.^2) + 0.5 * (control(terminal_time, control_params + uncert[i, :], num_params, terminal_time))^2 - gamma * sum(uncert[i, :].^2) <=  slack)

#     end
# end

function slack_constraint(model, terminal_state, lambda, gamma, control_params, slack, states_sampled, uncert, num_params, terminal_time)
    
    for i in 1:num_samples
    
        states_sampled_error = states_sampled[:, i] - terminal_state
        @constraint(model, lambda * sum(control_params.^2) + 0.5 * 100 * sum(states_sampled_error.^2) - gamma * sum(uncert[i, :].^2) <=  slack)

    end
end

function dynamics_constraints(model, initial_state, state_matrix, control_matrix, exp_comp_matrix, states_sampled, control_params, uncert, num_samples, num_params)

    initial_distribution = initial_state

    A = state_matrix
    B = control_matrix
    A_inv = state_matrix_inv
    
    for i in 1:num_samples

    psuedo_state = initial_distribution
    state = (exp(A * terminal_time)) * psuedo_state
    state  = state + control_integral(terminal_time, B, state_eigen_values, state_eigen_vectors, control_params, num_params)
    state = state + uncert_integral(terminal_time, B, state_eigen_values, state_eigen_vectors, uncert[i,:], num_params)

    @constraint(model, states_sampled[:, i] .== state)

    end
end