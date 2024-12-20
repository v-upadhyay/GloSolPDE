using Ipopt
using LinearAlgebra
using StaticArrays

function exponential_state_matrix(eigen_values, eigen_vectors, t)

    D = Diagonal(exp.(eigen_values * t))
    V = eigen_vectors
    return V * D * inv(V)

end

function integral_exponential_state_matrix(eigen_values, eigen_vectors, t_fin, t)

    D = Diagonal(exp.(eigen_values * (t_fin - t)) ./ (-eigen_values))
    V = eigen_vectors
    return V * D * inv(V)
    
end

function control_constraints(model, control, num_params, dT)
    
    for i in 1:num_params

        @constraint(model, -10 <= control[i] <= 10)

    end
end


function slack_constraint(model, control, slack, displacement, uncert, dT)
    
    lambda = 0.01
    gamma = 50.0
    # terminal_distribution = [10 * exp(-(x - 0.8)^2) for x in 1:discrete_points]
    terminal_distribution = [0 for x in 1:discrete_points]
    for i in 1:num_samples
    
        displacement_error = displacement[:, i] - terminal_distribution
        @constraint(model, lambda * sum(control.^2) + 0.5 * sum(displacement[:, i].^2) - gamma * sum(uncert[i,:,:].^2) <=  slack)

    end
    @constraint(model, slack <= 5)
end

function dynamics_constraints(model, state_matrix, control_matrix, state_matrix_inv, displacement, control, uncert, num_samples, dT, num_params)

    initial_distribution = [10 * (x/discrete_points)^2 * (1-(x/discrete_points))^3 for x in 1:discrete_points]

    A = state_matrix
    B = control_matrix
    A_inv = state_matrix_inv

    for i in 1:num_samples

    psuedo_state = initial_distribution
    state = (exp(A * num_params * dT)) * psuedo_state
    for j in 1:num_params
    
        state = state + (exp(A * ((num_params + 1 - j)*dT)) - exp(A * ((num_params - j)*dT))) * A_inv * (B * control[j] + uncert[i, :, j])
    
    end

    @constraint(model, displacement[:, i] .== state)

    end
end



function terminal_constraints(model, state_matrix, control_matrix, state_matrix_inv, displacement, control, uncert, num_samples, dT, num_params)

    initial_distribution = [1 * (x/discrete_points)^2 * (1-(x/discrete_points))^3 for x in 1:discrete_points]


    A = state_matrix
    B = control_matrix
    A_inv = state_matrix_inv

    for i in 1:num_samples
        
        for j in 1:discrete_points
        
            @constraint(model, displacement[j, i] <= 0.1)
        
        end
    end
end



