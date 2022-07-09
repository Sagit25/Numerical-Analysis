# LU decomposition by Gauss Elimination
# Confirms A = LU

using LinearAlgebra

# LU decomposition
function LU_decomposition(L, U, n)
    for k = 1:n-1
        for i = k+1:n
            L[i, k] = U[i, k] / U[k, k]
            U[i, k] = 0
            U[i, k+1:n] = U[i, k+1:n] - L[i, k] * U[k, k+1:n]
        end
    end
    return
end

# Check answer: A = LU
function check_ans(A, L, U, n)
    for i = 1:n
        for j = 1:n
            if abs(A[i, j] - dot(L[i, :], U[:, j])) > 1
                println("Your program give an error! : ", (i, j));
                exit()
            end
        end
    end
end

# Initialize A: 3*3 matrix
n = 3
A = [[2.0 -1 0]; [-1 2 -1]; [0 -1 2]]
L = zeros(n, n)
for i = 1 : n
    L[i, i] = 1
end
U = copy(A)

LU_decomposition(L, U, n)
check_ans(A, L, U, n)

# Print answer
println("L = ", L)
println("U = ", U)
println("LU = ", L * U)