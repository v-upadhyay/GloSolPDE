using LinearAlgebra

function compute_alpha_matrix(A)
    n = size(A, 1)
    λ = eigvals(A)  # Compute eigenvalues of A
    
    # Construct Vandermonde matrix for α coefficients
    V = [λ[i]^k for i in 1:n, k in 0:n-1]
    
    return inv(V)  # This is the matrix that gets multiplied to exp(λt) to obtain α
end

A = diagm([1, 2, 3])  # Example matrix
alpha_matrix = inv(compute_alpha_matrix(A))
println("Alpha matrix:")
println(alpha_matrix)