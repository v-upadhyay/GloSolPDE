using QuadGK, Serialization

# Compute Gram matrix Q
function compute_Q(num_params::Int, terminal_time::Float64, time_period::Float64)
    Q = zeros(num_params, num_params)

    # function basis(t, i)
    #     if i == 1
    #         return 1.0
    #     elseif i <= (num_params + 1) ÷ 2
    #         return sin(2π * (i - 1) * t / time_period)
    #     else
    #         return cos(2π * (i - (num_params + 1) ÷ 2) * t / time_period)
    #     end
    # end

    function basis(t, i)
        if i==1
            return 1
        elseif i <= (num_params+1)/2
            return sin(2*pi*(i-1)*t/time_period)
        else
            return cos(2 * pi * (i - (num_params+1) / 2) * t / time_period)
        end
    end

    for i in 1:num_params
        for j in i:num_params
            integrand(t) = basis(t, i) * basis(t, j)
            Q[i, j], _ = quadgk(integrand, 0, terminal_time)
            Q[j, i] = Q[i, j]  # symmetry
        end
    end

    return Q
end

# Parameters
terminal_time = 10.0
time_period = 10 * terminal_time
num_params = 20

# Compute and serialize
Q = compute_Q(num_params, terminal_time, time_period)

eigenvalue_Q = eigen(Q).values
println("Eigenvalues of Q: ", eigenvalue_Q)
serialize("Q_matrix.dat", Q)
