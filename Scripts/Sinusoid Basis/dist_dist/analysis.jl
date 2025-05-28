using JuMP
using Ipopt
using Plots
using LinearAlgebra
using StaticArrays
using Dates
using Serialization


include("system.jl")


##########################################################################
########################## PDE Parameters ################################
##########################################################################

boundary_args = SVector(1, 1, 0, 0)                                                                                 # dirichlet_0, dirichlet_1, neumann_0, neumann_1
discrete_points = 50                                                                                               # Number of discrete points in spatial domain
terminal_time = 5                                                                                                   # Terminal time
initial_distribution = [30 * ((x/100)^3) * (1 - x/100)^3 for x in 1:100]          # Initial distribution for the PDE



##########################################################################
########################## System Matrices ###############################
##########################################################################

state_matrix = state_params()    
control_matrix = control_params()
eigen_values, eigen_vectors = eigen(state_matrix)                                           # Eigen values and vectors of the state matrix
state_matrix_inv = inv(state_matrix)                                                        # Inverse of the state matrix


##########################################################################
########################## SIP Parameters ################################
##########################################################################

num_samples = 101
num_params = 100
overall_samples = num_samples * discrete_points * num_params
dT = terminal_time / num_params


##########################################################################
###################### Trajectory Computation ############################
##########################################################################

function state_information(initial_state, control, discrete_points, state_matrix, control_matrix, dT, time)
    interval = floor(time/dT)
    interval = Int(interval)
    left_time = time - interval * dT

    state = (exp(state_matrix * time)) * initial_state
    for i in 1:interval
        state = state + (exp(state_matrix * (time - (i - 1) * dT)) - exp(state_matrix * (time - i * dT))) * inv(state_matrix) * (control_matrix * control[i])
    end

    if interval < length(control)
        state = state + (exp(state_matrix * left_time) - Matrix{Float64}(I, discrete_points, discrete_points)) * inv(state_matrix) * (control_matrix * control[interval + 1])
    end

    return state
end


function state_trajectory(initial_state, control, discrete_points, state_matrix, control_matrix, dT, sample_times)
    N = length(sample_times)
    trajectory = zeros(discrete_points, N)

    for i in 1:N
        trajectory[:, i] = state_information(initial_state, control, discrete_points, state_matrix, control_matrix, dT, sample_times[i])
    end


    return trajectory
end


function feedback_trajectory(time_samples, control_data, dT)
    N = length(time_samples)
    control = zeros(N)

    for i in 1:N
        control[i] = control_data[Int(floor(time_samples[i]/dT)) + 1]
    end

    return control
end


##########################################################################
##################### Trajectory Visualization ###########################
##########################################################################

simulation_samples = 100
sample_time = [terminal_time*(i)/simulation_samples for i in 1:simulation_samples]

control = deserialize(open("Recent_Control_Trajectory.dat"))
trajectory = state_trajectory(initial_distribution, control, discrete_points, state_matrix, control_matrix, dT, sample_time)


##### Surface Plot #####


plot(trajectory, st=:surface, xlabel="Time", ylabel="Length", zlabel="Temperature", title="Heat Distribution Along the Rod")

##### All in one Plot #####
# plot(trajectory, xlabel="Length", ylabel="Temperature", title="Heat Distribution Along the Rod", legend=false)

# feedback = feedback_trajectory(sample_time, control, dT)
# plot(sample_time, feedback, xlabel="Time", ylabel="Control", title="Control Trajectory", legend=false)

# serialize("recent_traj.dat", trajectory)