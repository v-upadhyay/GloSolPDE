#=
############################################################################################
	Robust Optimal Control for Linear Parabolic PDEs
	<main.jl>
	© vaibhav.u@iitb.ac.in

	This piece of code is written for generating near optimal solutions 
   to Robust Optimal Control Problem for Linear Systems in continuous 
   time regime only with terminal constraints. Extensions to be written
   later for intermediate constraints.

	At the moment this file is accompanied by 
	1. main.jl containing the simulated anneling driven maximization over 
      the uncertainty set.
	2. internal_minimization.jl containing the minimization of cost function over the control
      variables for given samples of uncertainty.
   3. constraints.jl containing the constraints for the internal_minimization 
      optimization problem (control, slack & dynamics constraints).
   4. objective_function.jl containing the cost function for the internal_minimization.

	The code is semi-optimized; full-scale optimization of the code should be
	done later.
############################################################################################
=#



run(`clear`)


using JuMP
using Ipopt
using CairoMakie
using LinearAlgebra
using StaticArrays
using Dates
using Serialization
using Optim
using Random
using StaticArrays
using BlackBoxOptim

println("Reading the files...")

include("objective_function.jl")
include("constraints.jl")
include("internal_minimization.jl")
include("system.jl")
include("uncert_projection.jl")
# include("exponential_comp.jl")

println("Files read successfully")


function optimal_control(N)

##########################################################################
########################## OCP Parameters ################################
##########################################################################

terminal_time = 0.5                                                                                             # Terminal time
time_period = 10.0 * terminal_time                                                                                 # Time period for the control
terminal_state = [0.0 for x in 1:N]                                                                            # Initial distribution for the PDE
initial_state = [5 * (1 - 2*sin(1.5*pi*x/N)) for x in 1:N]                                                     # Terminal distribution for the PDE
lambda = 0.00001                                                                                                   # Weight for the control cost
gamma = 0.001                                                                                             # Weight for the uncertainty cost

##########################################################################
########################## System Matrices ###############################
##########################################################################


state_dim = size(initial_state)[1]                                         # Dimension of the state space
state_matrix = state_params(N)                                              # State matrix
control_matrix = control_params(N)                                          # Control matrix
control_bounds = [-20, 5]                                                   # Control bounds


##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_params = 20
num_samples = num_params + 1
disturbance_upper_bound = 0.01
disturbance_lower_bound = -0.01
time_lower_bound = 0.0
time_upper_bound = terminal_time
uncert_radius = (disturbance_upper_bound - disturbance_lower_bound)/2 

single_psuedo_uncert_row_lower_bound = vcat(fill(disturbance_lower_bound, num_params))
single_psuedo_uncert_row_upper_bound = vcat(fill(disturbance_upper_bound, num_params))

# Flatten these bounds for all num_samples rows of psuedo_uncert
flat_psuedo_uncert_lower_bounds = vcat([single_psuedo_uncert_row_lower_bound for _ in 1:num_samples]...)
flat_psuedo_uncert_upper_bounds = vcat([single_psuedo_uncert_row_upper_bound for _ in 1:num_samples]...)

# Bounds for time_samples
time_samples_lower_bounds = fill(time_lower_bound, num_samples)
time_samples_upper_bounds = fill(time_upper_bound, num_samples)

# Combine to form the final flat bound vectors
sip_lower_bound = vcat([0.0], flat_psuedo_uncert_lower_bounds, time_samples_lower_bounds)
sip_upper_bound = vcat([0.0], flat_psuedo_uncert_upper_bounds, time_samples_upper_bounds)

##########################################################################
########################## Maximization Cost #############################
##########################################################################

function maximization_cost(uncert_samples)

   # Decompose time and uncert from initial_guess
   time_samples = uncert_samples[end-num_samples+1:end]  # Extract the last num_samples elements as time samples

   uncert_array = uncert_samples[2:end-num_samples]   # Extract the rest as uncertainty samples
   uncert_0 = uncert_samples[1]  # Extract the first element as the uncertainty coefficient
   uncert_matrix = reshape(uncert_array, num_samples, num_params)  # Reshape into uncertainty tensor
   #  uncert = uncert_projector(uncert_matrix, uncert_radius, uncert_0, num_samples, num_params)  # Project the uncertainty samples onto the hyperplane 
   uncert = uncert_matrix

    cost = -internal_minimization(state_matrix, control_matrix, control_bounds, num_samples, num_params, terminal_time, uncert, state_dim, initial_state, terminal_state, gamma, lambda, time_samples, time_period)  # Minimize the cost function over the control variables
    
    return cost

end

##########################################################################
########################## Global Maximization ###########################
##########################################################################



# Run DE optimization
result = bboptimize(maximization_cost;
SearchRange = [ (lo, hi) for (lo, hi) in zip(sip_lower_bound, sip_upper_bound) ],
NumDimensions = length(sip_lower_bound),
Method = :de_rand_1_bin,
MaxSteps = 50
)


best_sol = best_candidate(result)
# best_val = best_fitness(result)

# Decompose time and uncert from initial_guess
opt_time_samples = best_sol[end-num_samples+1:end]  # Extract the last num_samples elements as time samples
opt_uncert_array = best_sol[2:end-num_samples]   # Extract the rest as uncertainty samples
contant_check = best_sol[1]  # Extract the first element as the uncertainty coefficient
println(contant_check)
# Reshape the uncertainty samples into the uncertainty tensor
opt_uncert = reshape(opt_uncert_array, num_samples, num_params)  # Reshape into uncertainty tensor
opt_uncert = uncert_projector(opt_uncert, uncert_radius, contant_check, num_samples, num_params)  # Project the uncertainty samples onto the hyperplane
serialize("disturbance_maximizer.dat", opt_uncert)
control_commands = control_trajectory(state_matrix, control_matrix, control_bounds, num_samples, num_params, terminal_time, opt_uncert, state_dim, initial_state, terminal_state, gamma, lambda, opt_time_samples, time_period)
serialize("control_trajectory.dat", control_commands)
println("New best solution saved for N:", N)

end



N_grid = 50

optimal_control(N_grid)
