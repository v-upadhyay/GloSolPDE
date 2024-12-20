





function internal_minimization(state_matrix, control_matrix, state_matrix_inv, num_samples, num_params, terminal_time, uncert_samples, discrete_points)

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
@variable(model, displacement[1:discrete_points, 1:num_samples])                   # Terminal distribution for each disturbance sample

# Define the objective function to minimize the total distance
cost_function(model, slack)

# Add constraints to the model
slack_constraint(model, control, slack, displacement, uncert_samples, dT)
dynamics_constraints(model, state_matrix, control_matrix, state_matrix_inv, displacement, control, uncert_samples, num_samples, dT, num_params)
control_constraints(model, control, num_params, dT)



# Solve the optimization problem
optimize!(model)

return objective_value(model)

end


function control_trajectory(state_matrix, control_matrix, state_matrix_inv, num_samples, num_params, terminal_time, uncert_samples, discrete_points)

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
    @variable(model, displacement[1:discrete_points, 1:num_samples])                   # Terminal distribution for each disturbance sample
    
    # Define the objective function to minimize the total distance
    cost_function(model, slack)
    
    # Add constraints to the model
    slack_constraint(model, control, slack, displacement, uncert_samples, dT)
    dynamics_constraints(model, state_matrix, control_matrix, state_matrix_inv, displacement, control, uncert_samples, num_samples, dT, num_params)
    control_constraints(model, control, num_params, dT)
    
    
    
    # Solve the optimization problem
    optimize!(model)
    
    return value.(control)
    
    end