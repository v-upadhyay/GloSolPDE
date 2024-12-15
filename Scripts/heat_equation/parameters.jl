
function system_matrix(discrete_points, boundary_args, parameters)

    theta = parameters[1]
    alpha = parameters[2]

    dirichlet_0 = boundary_args[1]
    dirichlet_1 = boundary_args[2]
    neumann_0 = boundary_args[3]
    neumann_1 = boundary_args[4]

    h = 1 / (discrete_points + 1)

    r_0 = neumann_0/ (3 * neumann_0  - 2 * h * dirichlet_0)
    r_1 = neumann_1/ (3 * neumann_1  + 2 * h * dirichlet_1)

    state_matrix = zeros(discrete_points, discrete_points)

    for i in 1:discrete_points
        state_matrix[i, i] = -2
        if i > 1
            state_matrix[i, i-1] = 1
        end
        if i < discrete_points
            state_matrix[i, i+1] = 1
        end
    end

    state_matrix[1, 1] = -2 + 4 * r_0
    state_matrix[1, 2] = 1 - r_0
    state_matrix[end, end] = -2 + 4 * r_1
    state_matrix[end, end-1] = 1 - r_1

    state_matrix = theta * state_matrix / h^2

    state_matrix = state_matrix + alpha * Matrix{Float64}(I, discrete_points, discrete_points)


    h = 1 / (discrete_points + 1)

    b_n = 2 * h * theta / (3 * neumann_1 + 2 * h * dirichlet_1)

    control_matrix = zeros(discrete_points, 1)

    control_matrix[end] = b_n / h^2

    return state_matrix, control_matrix

end
