#=
############################################################################################
	Robust Optimal Control for Linear Systems with Terminal Constraints
	<main.jl>
	© vaibhav.u@iitb.ac.in

	This piece of code is written for generating near optimal solutions 
   to Robust Optimal Control Problem for Linear Systems in continuous 
   time regime only with terminal constraints. Extensions to be written
   later for intermediate constraints.

	At the moment this file is accompanied by 
	1. global_maximization.jl containing the simulated anneling driven maximization over 
      the uncertainty set.
	2. internal_minimization.jl containing the minimization of cost function over the control
      variables for given samples of uncertainty.
   3. constraints.jl containing the constraints for the internal_minimization 
      optimization problem (control, slack & dynamics constraints).
   4. objective_function.jl containing the cost function for the internal_minimization.
   5. parameters.jl containing the system matrices for the linear system based on given 
      parameters and boundary conditions.

	This piece of code works in all moderately high dimensional systems: BELOW dim 100.
	The code is semi-optimized; full-scale optimization of the code should be
	done later.
############################################################################################
=#






using JuMP
using Ipopt
using Plots
using LinearAlgebra
using StaticArrays
using Dates
using Serialization
using Optim
using Random
using StaticArrays

println("Reading the files...")

include("objective_function.jl")
include("terminal_conditions.jl")
include("constraints.jl")
include("internal_minimization.jl")
include("parameters.jl")

println("Files read successfully")


##########################################################################
########################## PDE Parameters ################################
##########################################################################

boundary_args = SVector(1, 1, 0, 0)                                             # dirichlet_0, dirichlet_1, neumann_0, neumann_1
discrete_points = 100                                                           # Number of discrete points in spatial domain
parameters = SVector(1, 0.2 + pi ^ 2)                                           # theta, alpha
terminal_time = 1                                                               # Terminal time
initial_distribution = [10 * (x/discrete_points)^2 * (1-(x/discrete_points))^3 for x in 1:discrete_points]          # Initial distribution for the PDE



##########################################################################
########################## System Matrices ###############################
##########################################################################

state_matrix, control_matrix = system_matrix(discrete_points, boundary_args, parameters)     # Control matrix
eigen_values, eigen_vectors = eigen(state_matrix)                                           # Eigen values and vectors of the state matrix
state_matrix_inv = inv(state_matrix)                                                        # Inverse of the state matrix



##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_samples = 20
num_params = 400
overall_samples = num_samples * discrete_points * num_params
initial_guess = vcat([0.01 for i in 1:overall_samples])     # Initial guess for the SIP
sip_lower_bound = vcat([-0.01 for _ in 1:overall_samples])
sip_upper_bound = vcat([0.01 for _ in 1:overall_samples])

##########################################################################
########################## Global Maximization ###########################
##########################################################################

function maximization_cost(uncert_samples)

   uncert = reshape(uncert_samples, num_samples, discrete_points, num_params)   # Reshape the uncertainty samples into uncertainty tensor
   cost = -internal_minimization(state_matrix, control_matrix, state_matrix_inv, eigen_values, eigen_vectors, 
   num_samples, num_params, terminal_time, uncert, discrete_points)
   return cost

end

start_time = time()
max_anneling_iterations = 30
options = Optim.Options(iterations = max_anneling_iterations, store_trace = true, show_trace = true)
result = optimize(maximization_cost, sip_lower_bound, sip_upper_bound, initial_guess, SAMIN(nt = 20, ns = 10, rt = 0.9), options)

end_time = time()
println(end_time - start_time)

# Extract the trace of the optimization process
trace = result.trace

# Extract the objective values from the trace
objective_values = [trace[i].value for i in 1:length(trace)]

# Plot the objective values
plot(-objective_values, title="Simulated Annealing Optimization", xlabel="Iteration", ylabel="Objective Value", legend=false, ylim=(0, 0.02))