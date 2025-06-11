using LinearAlgebra, CairoMakie

include("system.jl")  # This should define state_params(N)

k = 80.2
D = 0.025
L = 0.2
A = pi * D^2 / 4

l = k * A / L

# Parameters
N = 10
T = 5
sampling_points = 50

# Matrix and eigen decomposition
A = state_params(N)
l = eigen(A).values
state_matrix_decomposition = eigen(A)
state_eigen_values = state_matrix_decomposition.values
state_eigen_vectors = state_matrix_decomposition.vectors
# println("Max eigenvalue: ", maximum(l))\
# println("Eigenvalues: ", diagm(state_eigen_values))

A_1 = inv(state_eigen_vectors) * exp(5*Diagonal(state_eigen_values)) * (state_eigen_vectors)

# println("Error: ", maximum(A_1-exp(5*A)))

l = zeros(10)
for i in 1:10
    l[i] = 2*rand() - 1
end
println("Random vector: ", abs.(l))
nor_l = l / sum(abs.(l))
println("Normalized vector: ", nor_l)
println("Sum: ", sum(abs.(nor_l)))


# # Initial condition and space-time grid
# x0 = [10 * sin(i * π / N) for i in 1:N]
# space = LinRange(0, 1, N)
# time = LinRange(0, T, sampling_points + 1)

# # Simulate evolution X(t) = exp(At) x₀
# X = zeros(length(time), N)
# for (i, tᵢ) in enumerate(time)
#     X[i, :] = exp(A * tᵢ) * x0
# end
# println("Check1")
# # Set up CairoMakie figure
# fig = Figure(size = (900, 640), backgroundcolor = :white)

# label_size = 28
# tick_size = 16
# println("Check2")
# ax = Axis3(fig[1, 1],
#     xlabel = "Space", ylabel = "Time", zlabel = "Amplitude",
#     xlabelsize = label_size, ylabelsize = label_size, zlabelsize = label_size,
#     xticklabelsize = tick_size, yticklabelsize = tick_size, zticklabelsize = tick_size,
#     xlabelcolor = :black, ylabelcolor = :black, zlabelcolor = :black,
#     xgridvisible = true, ygridvisible = true, zgridvisible = true,
#     xgridcolor = (:gray, 0.3), ygridcolor = (:gray, 0.3), zgridcolor = (:gray, 0.3),
#     xgridwidth = 1, ygridwidth = 1, zgridwidth = 1,
#     backgroundcolor = :white,
#     title = "Evolution of State over Time",
#     titlesize = 20
# )
# println("Check3")
# # Plot surface
# surface!(ax, space, time, X', colormap = :inferno, shading = true)
# println("Check4")
# # Save and show
# save("surface_plot.png", fig)
# fig