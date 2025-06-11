





function internal_minimization(state_matrix, control_matrix, control_bounds, num_samples, num_params, terminal_time, uncert, state_dim, initial_state, terminal_state, gamma, lambda, time_samples, time_period)


# Define the optimization model
model = Model(Ipopt.Optimizer)

# Set the maximum number of iterations
set_optimizer_attribute(model, "print_level", 0) 
set_optimizer_attribute(model, "max_iter", 400)
set_optimizer_attribute(model, "tol", 1e-3)
set_optimizer_attribute(model, "acceptable_tol", 1e-3)
set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
set_optimizer_attribute(model, "mu_strategy", "adaptive")
set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
                       

########################################
######### Control Variables ############
########################################
@variable(model, control_params[1:num_params])                                            # piecewise constant control
@variable(model, slack)                                                            # slack variable
@variable(model, states_sampled[1:state_dim, 1:num_samples])                   # Terminal distribution for each disturbance sample

# Define the objective function to minimize the total distance
cost_function(model, slack)

# Add constraints to the model
slack_constraint(model, terminal_state, lambda, gamma, control_params, slack, states_sampled, uncert, num_params, terminal_time)
dynamics_constraints(model, initial_state, state_matrix, control_matrix, states_sampled, control_params, uncert, num_samples, num_params, time_period, terminal_time)
control_constraints(model, control_params, uncert, num_params, control_bounds, time_samples, time_period)


# Solve the optimization problem
optimize!(model)

error_array = value.(states_sampled)
error = [norm(error_array[:, i]) for i in 1:num_samples]
println("Maximum Error Norm: ", maximum(error))
println("Solution status:" , termination_status(model))
return objective_value(model)

end


function control_trajectory(state_matrix, control_matrix, control_bounds, num_samples, num_params, terminal_time, uncert, state_dim, initial_state, terminal_state, gamma, lambda, time_samples, time_period)


    # Define the optimization model
    model = Model(Ipopt.Optimizer)
    
    # Set the maximum number of iterations
    # set_optimizer_attribute(model, "print_level", 0) 
    set_optimizer_attribute(model, "max_iter", 400)
    set_optimizer_attribute(model, "tol", 1e-3)
    set_optimizer_attribute(model, "acceptable_tol", 1e-3)
    set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    set_optimizer_attribute(model, "mu_strategy", "adaptive")
    set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
                           
    
    ########################################
    ######### Control Variables ############
    ########################################
    @variable(model, control_params[1:num_params])                                            # piecewise constant control
    @variable(model, slack)                                                            # slack variable
    @variable(model, states_sampled[1:state_dim, 1:num_samples])                   # Terminal distribution for each disturbance sample
    
    # Define the objective function to minimize the total distance
    cost_function(model, slack)
    
    # Add constraints to the model
    slack_constraint(model, terminal_state, lambda, gamma, control_params, slack, states_sampled, uncert, num_params, terminal_time)
    dynamics_constraints(model, initial_state, state_matrix, control_matrix, states_sampled, control_params, uncert, num_samples, num_params, time_period, terminal_time)
    control_constraints(model, control_params, uncert, num_params, control_bounds, time_samples, time_period)
    
    
    # Solve the optimization problem
    optimize!(model)

    error_array = value.(states_sampled)
    error = [norm(error_array[:, i]) for i in 1:num_samples]
    println("Maximum Error Norm: ", maximum(error))
    
    return value.(control_params)
    
    end