using LinearAlgebra

function theta(x)

    return 1

end

function sigma(x)

    return 2

end

function alpha(x)

    return 0.3

end



function state_params(discrete_points)

    h = 1 / (discrete_points + 1)

    theta_n = diagm([theta(x * h) for x in 1:discrete_points])
    sigma_n = diagm([sigma(x * h) for x in 1:discrete_points])
    alpha_n = diagm([alpha(x * h) for x in 1:discrete_points])

    dirichlet_0 = 1
    dirichlet_1 = 0
    neumann_0 = sigma(0)/2
    neumann_1 = 1

    r_0 = neumann_0/ (3 * neumann_0  - 2 * h * dirichlet_0)
    r_1 = neumann_1/ (3 * neumann_1  + 2 * h * dirichlet_1)
    q_0 = -dirichlet_0 / (neumann_0 - h * dirichlet_0)

    state_matrix = zeros(discrete_points, discrete_points)

    L_n = zeros(discrete_points, discrete_points)

    for i in 1:discrete_points
        L_n[i, i] = -2
        if i > 1
            L_n[i, i-1] = 1
        end
        if i < discrete_points
            L_n[i, i+1] = 1
        end
    end

    L_n[1, 1] = -2 + 4 * r_0
    L_n[1, 2] = 1 - r_0
    L_n[end, end] = -2 + 4 * r_1
    L_n[end, end-1] = 1 - r_1

    L_n = L_n / h^2

    D_n = zeros(discrete_points, discrete_points)

    for i in 2:discrete_points
        D_n[i, i-1] = -1
        D_n[i, i] = 1
    end

    D_n[1, 1] = h * q_0

    D_n = D_n / h

    state_matrix = theta_n * L_n + sigma_n * D_n + alpha_n

    return state_matrix

end


function control_params(discrete_points)

    dirichlet_0 = 1
    dirichlet_1 = 0
    neumann_0 = sigma(0)/2
    neumann_1 = 1

    h = 1 / (discrete_points + 1)

    b_n = 2 * h * theta(1) / (3 * neumann_1 + 2 * h * dirichlet_1)

    control_matrix = zeros(discrete_points, 1)

    control_matrix[end] = b_n / h^2

    return control_matrix

end



