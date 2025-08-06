using LinearAlgebra

function theta(x)

    k_steel = 0.188

    return k_steel

end

function sigma(x)

    return 0

end

function alpha(x)

    return 0

end



function state_params(discrete_points)

    h = 1 / (discrete_points + 1)

    theta_n = diagm([theta(x * h) for x in 1:discrete_points])
    sigma_n = diagm([sigma(x * h) for x in 1:discrete_points])
    alpha_n = diagm([alpha(x * h) for x in 1:discrete_points])

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

    L_n[1, 1] = -5/2
    L_n[1, 2] = 2
    L_n[1, 3] = -1/2
    L_n[end, end] = -5/2
    L_n[end, end-1] = 2
    L_n[end, end-2] = -1/2

    L_n = L_n / h^2

    state_matrix = theta_n * L_n

    return state_matrix

end


function control_params(discrete_points)

    h = 1 / (discrete_points + 1)

    control_matrix = zeros(discrete_points, 1)

    control_matrix[end] = theta(1) / h^2

    return control_matrix

end



