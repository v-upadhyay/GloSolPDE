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
	1. main.jl containing the simulated anneling driven maximization over 
      the uncertainty set.
	2. internal_minimization.jl containing the minimization of cost function over the control
      variables for given samples of uncertainty.
   3. constraints.jl containing the constraints for the internal_minimization 
      optimization problem (control, slack & dynamics constraints).
   4. objective_function.jl containing the cost function for the internal_minimization.
   5. parameters.jl containing the system matrices for the linear system based on given 
      parameters and boundary conditions.
   6. analysis.jl containing the verification tools for the validity of the solution through
      trajectory visualization in the nominal case (zeros disturbance).

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

boundary_args = SVector(1, 1, 0, 0)                                                                                 # dirichlet_0, dirichlet_1, neumann_0, neumann_1
discrete_points = 100                                                                                               # Number of discrete points in spatial domain
parameters = SVector(1, 0.2 + pi ^ 2)                                                                               # theta, alpha
terminal_time = 1                                                                                                   # Terminal time
initial_distribution = [10 * (x/discrete_points)^2 * (1-(x/discrete_points))^3 for x in 1:discrete_points]          # Initial distribution for the PDE



##########################################################################
########################## System Matrices ###############################
##########################################################################

state_matrix, control_matrix = system_matrix(discrete_points, boundary_args, parameters)    # Control matrix
eigen_values, eigen_vectors = eigen(state_matrix)                                           # Eigen values and vectors of the state matrix
state_matrix_inv = inv(state_matrix)                                                        # Inverse of the state matrix



##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_samples = 51
num_params = 50
overall_samples = num_samples * discrete_points * num_params
initial_guess = vcat([0.1 * (-1)^i for i in 1:overall_samples])     # Initial guess for the SIP
# initial_guess = deserialize(open("disturbance.dat"))
sip_lower_bound = vcat([-0.1 for _ in 1:overall_samples])
sip_upper_bound = vcat([0.1 for _ in 1:overall_samples])

##########################################################################
########################## Maximization Cost #############################
##########################################################################

function maximization_cost(uncert_samples)

   uncert = reshape(uncert_samples, num_samples, discrete_points, num_params)   # Reshape the uncertainty samples into uncertainty tensor
   cost = -internal_minimization(state_matrix, control_matrix, state_matrix_inv, num_samples, num_params, terminal_time, uncert, discrete_points)
   return cost

end


##########################################################################
########################## Global Maximization ###########################
##########################################################################

anneling_trials = 20                   # Number of anneling trials from same initial guess                                                      
max_anneling_iterations = 40           # Number of iterations anneling in each trial




options = Optim.Options(iterations = max_anneling_iterations, store_trace = true, show_trace = true)
result = optimize(maximization_cost, sip_lower_bound, sip_upper_bound, initial_guess, SAMIN(nt = 10, ns = 100, rt = 0.95), options)

objective = zeros(anneling_trials + 1, max_anneling_iterations)

# Extract the trace of the optimization process
trace = result.trace
optimal_disturbnace = result.minimizer

# Extract the objective values from the trace
objective[1, :] = [trace[i].value for i in 1:length(trace)]


for i in 1:anneling_trials

start_time = time()
options = Optim.Options(iterations = max_anneling_iterations, store_trace = true, show_trace = true)
result = optimize(maximization_cost, sip_lower_bound, sip_upper_bound, initial_guess, SAMIN(nt = 10, ns = 100, rt = 0.95), options)

end_time = time()
println(end_time - start_time)

# Extract the trace of the optimization process
trace = result.trace
optimal_disturbnace = result.minimizer

# Extract the objective values from the trace
objective[i+1, :] = [trace[i].value for i in 1:length(trace)]

end

serialize("multi_objective_runs.dat", objective)
serialize("disturbance_maximizer.dat", optimal_disturbnace)

control_commands = control_trajectory(state_matrix, control_matrix, state_matrix_inv, num_samples, num_params, terminal_time, reshape(optimal_disturbnace, num_samples, discrete_points, num_params), discrete_points)
serialize("Recent_Control_Trajectory.dat", control_commands)