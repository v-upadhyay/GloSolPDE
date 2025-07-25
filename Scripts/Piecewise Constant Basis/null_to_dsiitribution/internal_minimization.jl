





function internal_minimization(state_matrix, control_matrix, state_matrix_inv, control_bounds, num_samples, num_params, terminal_time, uncert_samples, state_dim, initial_state, terminal_state, gamma, lambda)

dT = terminal_time / num_params

# Define the optimization model
model = Model(Ipopt.Optimizer)

# Set the maximum number of iterations
set_optimizer_attribute(model, "print_level", 0) 
set_optimizer_attribute(model, "max_iter", 100)
set_optimizer_attribute(model, "tol", 1e-3)
set_optimizer_attribute(model, "acceptable_tol", 1e-3)
set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
set_optimizer_attribute(model, "mu_strategy", "adaptive")
set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
                       

########################################
######### Control Variables ############
########################################
@variable(model, control[1:num_params])                                            # piecewise constant control
@variable(model, slack)                                                            # slack variable
@variable(model, states_sampled[1:state_dim, 1:num_samples])                   # Terminal distribution for each disturbance sample

# Define the objective function to minimize the total distance
cost_function(model, slack)

# Add constraints to the model
slack_constraint(model, terminal_state, lambda, gamma, control, slack, states_sampled, uncert_samples, dT)
dynamics_constraints(model, initial_state, state_matrix, control_matrix, state_matrix_inv, states_sampled, control, uncert_samples, num_samples, dT, num_params)
control_constraints(model, control, num_params, control_bounds, dT)



# Solve the optimization problem
optimize!(model)


terminal_array = value.(states_sampled)
error_array = zeros(state_dim, num_samples)
for i in 1:num_samples
    error_array[:, i] = terminal_array[:, i] - terminal_state
end
error = [norm(error_array[:, i]) for i in 1:num_samples]
print("Maximum Error Norm: ", maximum(error))

return objective_value(model)

end


function control_trajectory(state_matrix, control_matrix, state_matrix_inv, control_bounds, num_samples, num_params, terminal_time, uncert_samples, state_dim, initial_state, terminal_state, gamma, lambda)

    dT = terminal_time / num_params
    
    # Define the optimization model
    model = Model(Ipopt.Optimizer)
    
    # Set the maximum number of iterations
    set_optimizer_attribute(model, "print_level", 0) 
    set_optimizer_attribute(model, "max_iter", 100)
    set_optimizer_attribute(model, "tol", 1e-3)
    set_optimizer_attribute(model, "acceptable_tol", 1e-3)
    set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    set_optimizer_attribute(model, "mu_strategy", "adaptive")
    set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
                           
    
    ########################################
    ######### Control Variables ############
    ########################################
    @variable(model, control[1:num_params])                                            # piecewise constant control
    @variable(model, slack)                                                            # slack variable
    @variable(model, states_sampled[1:state_dim, 1:num_samples])                   # Terminal distribution for each disturbance sample
    
    # Define the objective function to minimize the total distance
    cost_function(model, slack)
    
    # Add constraints to the model
    slack_constraint(model, terminal_state, lambda, gamma, control, slack, states_sampled, uncert_samples, dT)
    dynamics_constraints(model, initial_state, state_matrix, control_matrix, state_matrix_inv, states_sampled, control, uncert_samples, num_samples, dT, num_params)
    control_constraints(model, control, num_params, control_bounds, dT)
    
    
    
    # Solve the optimization problem
    optimize!(model)

    terminal_condition = value.(states_sampled)
    plot(terminal_condition)
    plot!(terminal_state)
    
    return value.(control)
    
    end