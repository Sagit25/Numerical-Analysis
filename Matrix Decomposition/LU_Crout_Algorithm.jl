# LU decomposition by Crout algorithm
# Solve Ax = b for Tridiagonal Matrix

using LinearAlgebra

# Crout algorithm
function crout_algorithm(A, L, U, n)
    for i = 1:n
        for j = i:n
            L[j, i] = A[j, i] - dot(L[j, 1:i-1], U[1:i-1, i])
        end
        for j = i+1:n
            U[i, j] = A[i, j] - dot(L[i, 1:i-1], U[1:i-1, j])
            U[i, j] = U[i, j] / L[i, i]
        end
    end
    return
end

# Forward elimination for lower triangular
function forward_elimination_lower(L, b, n)
    for j = 1:n
        b[j] = (b[j] - dot(L[j, 1:j-1], b[1:j-1])) / L[j, j]
    end
    return
end

# Back substitution for unit upper triangular
function back_substitution_unit_upper(U, b, n)
    for j = n:-1:1
        b[j] = b[j] - dot(U[j, j+1:n], b[j+1:n]);
    end
    return
end

# Check answer to calculate e = b-Ax
function check_ans(savea, saveb, x, n, err)
    for j = 1:n
        if abs(saveb[j] - dot(savea[j, 1:n], x[1:n])) > err
            println("Solution is wrong!")
            exit()
        end
    end
    return
end

# Initialize A, b: Using tridiagonal matrix 
print("Dimension: n = ")
n = parse(UInt64, readline())
err = 5.e-13
A = zeros(n, n)
for i = 1:n
    A[i, i] = 2
end
for j = 1:n-1
    A[j, j+1] = A[j+1, j] = -1
end
b = zeros(n)
b[1] = 1
L = zeros(n, n) 
U = zeros(n, n)
for i = 1:n
    U[i, i] = 1
end

# Copy A, b for check
savea = copy(A) 
saveb = copy(b)

crout_algorithm(A, L, U, n)
forward_elimination_lower(L, b, n)
back_substitution_unit_upper(U, b, n)
check_ans(savea, saveb, b, n, err)

println("Solution: x = " )
println(b)