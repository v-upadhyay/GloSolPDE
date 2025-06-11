using LinearAlgebra
using Ipopt


function uncert_projector(uncert_matrix, radius, c, num_samples, num_params)

    b = zeros(num_params)
    for i in 1:num_params
        if i == 1
            b[i] = 0
        elseif 1 < i <= (num_params+1)/2
            b[i] = 0
        else
            b[i] = 0
        end
    end

    solution_matrix = zeros(size(uncert_matrix))
    num_params = length(b)


    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0) 
    set_optimizer_attribute(model, "max_iter", 100) 
    set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
    # set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    @variable(model, x[1:num_params])

    @constraint(model, -radius .<= x[1:num_params] .<= radius)

    @constraint(model, dot(x, b) == c)


    for i in 1:num_samples
        @objective(model, Min, sum((x-uncert_matrix[i, :]).^2))
        optimize!(model)

        solution_matrix[i, :] = value.(x)
    end

    return solution_matrix
end


# num_params = 20



# c = 0.0
# num_samples = 21
# uncert_array_test = zeros(num_params, num_samples)
# for i in 1:num_samples
#     uncert_array_test[:, i] = [2*(rand() - 0.5) for i in 1:num_params]
# end
# radius = 1.0

# elapsed_time = @elapsed solution = uncert_projector(uncert_array_test, radius, c, num_samples, num_params)
# println("Time taken: $elapsed_time seconds")

# println("Solution:" , dot(b, solution[:, i]))