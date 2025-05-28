using Ipopt
using LinearAlgebra
using StaticArrays


function control_constraints(model, control, num_params, control_bounds, dT)
    
    lower_bound = control_bounds[1]
    upper_bound = control_bounds[2]
    epsilon = 0.001

    for i in 1:num_params

        @constraint(model, lower_bound <= control[i] <= upper_bound)

    end

    for i in 1:num_params-1

        @constraint(model, -epsilon <= control[i+1] - control[i] <= epsilon)

    end
end


function slack_constraint(model, terminal_state, lambda, gamma, control, slack, states_sampled, uncert, dT)
    
    for i in 1:num_samples
    
        states_sampled_error = states_sampled[:, i] - terminal_state
        @constraint(model, 0.5 * (lambda * sum(control.^2) * dT + sum(states_sampled_error.^2) - gamma * sum(uncert[i,:,:].^2) * dT) <=  slack)

    end
end

function dynamics_constraints(model, initial_state, state_matrix, control_matrix, state_matrix_inv, states_sampled, control, uncert, num_samples, dT, num_params)

    initial_distribution = initial_state

    A = state_matrix
    B = control_matrix
    A_inv = state_matrix_inv
    
    for i in 1:num_samples

    psuedo_state = initial_distribution
    state = (exp(A * num_params * dT)) * psuedo_state
    for j in 1:num_params
    
        state = state + (exp(A * ((num_params + 1 - j)*dT)) - exp(A * ((num_params - j)*dT))) * A_inv * (B * control[j] + uncert[i, :, j])
    
    end

    @constraint(model, states_sampled[:, i] .== state)

    end
end



