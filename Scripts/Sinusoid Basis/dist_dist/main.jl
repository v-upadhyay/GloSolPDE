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

	This piece of code works in all moderately high dimensional systems: BELOW dim 90.
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
include("exponential_comp.jl")

println("Files read successfully")


##########################################################################
########################## OCP Parameters ################################
##########################################################################

initial_modes = 3
terminal_mode = 1
terminal_time = 2.0                                                                                             # Terminal time
N = 50                                                                                                         # Semi-discretization grid points
terminal_state = [4 * ((x/N)) * (1 - x/N) for x in 1:N]                                                                            # Initial distribution for the PDE
initial_state = initial_state = [100 * ((x/N)^2) * (1 - x/N)^3 * sin(initial_modes * pi * x/N) for x in 1:N]                                                     # Terminal distribution for the PDE
lambda = 0.1                                                                                                   # Weight for the control cost
gamma = 0.5                                                                                                  # Weight for the uncertainty cost

##########################################################################
########################## System Matrices ###############################
##########################################################################

state_dim = size(initial_state)[1]                                         # Dimension of the state space
state_matrix = state_params(N)                                              # State matrix
state_matrix_decomposition = eigen(state_matrix)
state_eigen_values = state_matrix_decomposition.values
state_eigen_vectors = state_matrix_decomposition.vectors
exp_comp_matrix = compute_alpha_matrix(state_matrix)
control_matrix = control_params(N)                                          # Control matrix
state_matrix_inv = inv(state_matrix)                                       # Inverse of the state matrix
control_bounds = [-5, 5]                                                   # Control bounds


##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_params = 20
num_samples = num_params + 1
anneling_trials = 0                   # Number of annealing trials from the same initial guess                                                      
max_anneling_iterations = 50          # Number of iterations annealing in each trial
overall_samples = num_samples * num_params
disturbance_upper_bound = 0.1
disturbance_lower_bound = -0.1
time_lower_bound = 0.0
time_upper_bound = terminal_time

single_psuedo_uncert_row_lower_bound = vcat([disturbance_lower_bound], fill(-1.0, num_params))
single_psuedo_uncert_row_upper_bound = vcat([disturbance_upper_bound], fill(1.0, num_params))

# Flatten these bounds for all num_samples rows of psuedo_uncert
flat_psuedo_uncert_lower_bounds = vcat([single_psuedo_uncert_row_lower_bound for _ in 1:num_samples]...)
flat_psuedo_uncert_upper_bounds = vcat([single_psuedo_uncert_row_upper_bound for _ in 1:num_samples]...)

# Bounds for time_samples
time_samples_lower_bounds = fill(time_lower_bound, num_samples)
time_samples_upper_bounds = fill(time_upper_bound, num_samples)

# Combine to form the final flat bound vectors
sip_lower_bound = vcat(flat_psuedo_uncert_lower_bounds, time_samples_lower_bounds)
sip_upper_bound = vcat(flat_psuedo_uncert_upper_bounds, time_samples_upper_bounds)

# Lower and upper bounds for disturbance and time
# sip_lower_bound = vcat([disturbance_lower_bound for _ in 1:overall_samples], [time_lower_bound for _ in 1:num_samples])
# sip_upper_bound = vcat([disturbance_upper_bound for _ in 1:overall_samples], [time_upper_bound for _ in 1:num_samples])

##########################################################################
########################## Maximization Cost #############################
##########################################################################

function uncert_converter(uncert_array)

      # Reshape the uncertainty samples into the uncertainty tensor
      psuedo_uncert = reshape(uncert_array, num_samples, num_params+1)  # Reshape into uncertainty tensor
      uncert_coeff = psuedo_uncert[:, 1]  # Extract the first column as uncertainty coefficients
      uncert = zeros(num_samples, num_params)  # Initialize the uncertainty tensor
      for i in 1:num_samples
         denominator = sum(abs.(psuedo_uncert[i, 2:end]))  # Compute the denominator for each sample
         uncert[i, :] = uncert_coeff[i] * (psuedo_uncert[i, 2:end] ./ denominator)  # Normalize the uncertainty coefficients
         println("Uncertainty sum check $i: ", sum(abs.(uncert[i, :])), "\n")
      end
      return uncert
end

function maximization_cost(uncert_samples)

   # Decompose time and uncert from initial_guess
   time_samples = uncert_samples[end-num_samples+1:end]  # Extract the last num_samples elements as time samples
   uncert_array = uncert_samples[1:end-num_samples]   # Extract the rest as uncertainty samples

   uncert = uncert_converter(uncert_array)  # Convert the uncertainty samples into a tensor

   cost = -internal_minimization(state_matrix, control_matrix, exp_comp_matrix, control_bounds, num_samples, num_params, terminal_time, uncert, state_dim, initial_state, terminal_state, gamma, lambda, time_samples)  # Minimize the cost function over the control variables
   return cost

end

##########################################################################
########################## Global Maximization ###########################
##########################################################################

de_iterations = 1

current_best = 0.0

for i in 1:de_iterations

   # Run DE optimization
   result = bboptimize(maximization_cost;
   SearchRange = [ (lo, hi) for (lo, hi) in zip(sip_lower_bound, sip_upper_bound) ],
   NumDimensions = length(sip_lower_bound),
   Method = :de_rand_1_bin,
   MaxSteps = 2
   )


   best_sol = best_candidate(result)
   best_val = best_fitness(result)

   # if current_best <= - best_val
      global current_best = -best_val

      # Decompose time and uncert from initial_guess
      opt_time_samples = best_sol[end-num_samples+1:end]  # Extract the last num_samples elements as time samples
      opt_uncert_array = best_sol[1:end-num_samples]   # Extract the rest as uncertainty samples

      # Reshape the uncertainty samples into the uncertainty tensor
      opt_uncert = uncert_converter(opt_uncert_array)  # Convert the uncertainty samples into a tensor

      serialize("disturbance_maximizer.dat", opt_uncert)
      control_commands = control_trajectory(state_matrix, control_matrix, exp_comp_matrix, control_bounds, num_samples, num_params, terminal_time, opt_uncert, state_dim, initial_state, terminal_state, gamma, lambda, opt_time_samples)
      serialize("Recent_Control_Trajectory.dat", control_commands)
      println("New best solution saved")
   # end

end
