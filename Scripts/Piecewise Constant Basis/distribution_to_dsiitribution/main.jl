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

	This piece of code works in all moderately high dimensional systems: BELOW dim 100.
	The code is semi-optimized; full-scale optimization of the code should be
	done later.
############################################################################################
=#



run(`clear`)


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
include("constraints.jl")
include("internal_minimization.jl")
include("system.jl")

println("Files read successfully")


##########################################################################
########################## OCP Parameters ################################
##########################################################################

terminal_time = 10                                                                                             # Terminal time
terminal_state = [10 * ((x/100)^2) * (1 - x/100)^3 for x in 1:100]                                                                          # Initial distribution for the PDE
initial_state = [0.4 * sin(0.8 * pi * x/100) for x in 1:100]                                              # Terminal distribution for the PDE
lambda = 0.1                                                                                                   # Weight for the control cost
gamma = 500.0                                                                                                    # Weight for the uncertainty cost

##########################################################################
########################## System Matrices ###############################
##########################################################################

state_dim = size(initial_state)[1]                                         # Dimension of the state space
state_matrix = state_params()                                              # State matrix
control_matrix = control_params()                                          # Control matrix
state_matrix_inv = inv(state_matrix)                                       # Inverse of the state matrix
control_bounds = [-1, 1]                                                   # Control bounds


##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_params = 100
num_samples = num_params + 1
anneling_trials = 7                   # Number of anneling trials from same initial guess                                                      
max_anneling_iterations = 45           # Number of iterations anneling in each trial
overall_samples = num_samples * state_dim * num_params
disturbance_upper_bound = 0.01
disturbance_lower_bound = -0.01
initial_guess = vcat([disturbance_upper_bound * (-1)^i for i in 1:overall_samples])     # Initial guess for the SIP
sip_lower_bound = vcat([disturbance_lower_bound for _ in 1:overall_samples])
sip_upper_bound = vcat([disturbance_upper_bound for _ in 1:overall_samples])

##########################################################################
########################## Maximization Cost #############################
##########################################################################

function maximization_cost(uncert_samples)

   uncert = reshape(uncert_samples, num_samples, state_dim, num_params)   # Reshape the uncertainty samples into uncertainty tensor
   cost = -internal_minimization(state_matrix, control_matrix, state_matrix_inv, control_bounds, num_samples, num_params, terminal_time, uncert, state_dim, initial_state, terminal_state, gamma, lambda)  # Minimize the cost function over the control variables
   return cost

end


##########################################################################
########################## Global Maximization ###########################
##########################################################################

objective = zeros(anneling_trials + 1, max_anneling_iterations)

options = Optim.Options(iterations = max_anneling_iterations, store_trace = true, show_trace = true)
result = optimize(maximization_cost, sip_lower_bound, sip_upper_bound, initial_guess, SAMIN(nt = 10, ns = 100, rt = 0.8), options)

# Extract the trace of the optimization process
trace = result.trace
optimal_disturbnace = result.minimizer

# Extract the objective values from the trace
objective[1, :] = [trace[i].value for i in 1:length(trace)]


for i in 1:anneling_trials

   print("Ongoing Interation:", i)
   result = optimize(maximization_cost, sip_lower_bound, sip_upper_bound, initial_guess, SAMIN(nt = 10, ns = 100, rt = 0.95), options)

   # Extract the trace of the optimization process
   trace = result.trace
   optimal_disturbnace = result.minimizer

   # Extract the objective values from the trace
   objective[i+1, :] = [trace[i].value for i in 1:length(trace)]

end


serialize("objective_iteration_null_D2D.dat", objective)
serialize("disturbance_maximizer.dat", optimal_disturbnace)
control_commands = control_trajectory(state_matrix, control_matrix, state_matrix_inv, control_bounds, num_samples, num_params, terminal_time, reshape(optimal_disturbnace, num_samples, state_dim, num_params), state_dim, initial_state, terminal_state, gamma, lambda)
serialize("Recent_Control_Trajectory_D2D.dat", control_commands)
# plot(control_commands, title = "Control Trajectory", xlabel = "Time", ylabel = "Control", label = "Control", lw = 2, legend = :bottomright)
plot(-objective', title = "Objective Value vs Iteration", xlabel = "Iteration", ylabel = "Objective Value", label = "Objective Value", lw = 2, legend = :bottomright)