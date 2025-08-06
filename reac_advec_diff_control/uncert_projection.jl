using LinearAlgebra
using Ipopt


function basis(t, i, num_params, time_period)
    if i==1
        return 1
    elseif i <= (num_params+1)/2
        return sin(2*pi*(i-1)*t/time_period)
    else
        return cos(2*pi*(i-(num_params+1)/2)*t/time_period)
    end
end

function control(t, control_params, num_params, time_period)
    control_commands = 0.0
    for i in 1:num_params
        control_commands = control_commands + basis(t, i, num_params, time_period) * control_params[i]
    end
    return control_commands
end

function uncert_projector(uncert_matrix, radius, c, num_samples, num_params)

    solution_matrix = zeros(size(uncert_matrix))
    # num_params = length(b)

    for i in 1:num_samples

        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0) 
        set_optimizer_attribute(model, "max_iter", 100) 
        set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
        # set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
        @variable(model, x[1:num_params])

        @constraint(model, -radius .<= x[1:num_params] .<= radius)

        @constraint(model, control(0, x, num_params, 1) == c)

        @objective(model, Min, sum((x-uncert_matrix[i, :]).^2))
        optimize!(model)

        solution_matrix[i, :] = value.(x)
    end

    return solution_matrix
end


# # num_params = 20



# c = 0.0
# num_samples = 21
# uncert_array_test = zeros(num_samples, num_params)
# for i in 1:num_samples
#     uncert_array_test[i, :] = [2*(rand() - 0.5) for i in 1:num_params]
# end
# radius = 0.003

# elapsed_time = @elapsed solution = uncert_projector(uncert_array_test, radius, c, num_samples, num_params)
# println("Time taken: $elapsed_time seconds")
# error = [dot(b, solution[i, :]) for i in 1:num_samples]
# println("Solution error:" , maximum(error))