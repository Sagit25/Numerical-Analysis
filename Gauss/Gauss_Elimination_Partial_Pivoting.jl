# Gauss Elimination with partial pivoting
# Solve Ax = b for Matrix Need Partial Pivoting

using LinearAlgebra

# Partial pivoting
function partial_pivoting(a, b, n, j)
    js = argmax(abs.(a[j:n, j]))+j-1
    tmp = copy(a[j, :])
    a[j, :] = copy(a[js, :])
    a[js, :] = copy(tmp)
    tmp[j] = copy(b[j])
    b[j] = copy(b[js])
    b[js] = copy(tmp[j])
    return
end

# Gauss-Elimination with partial pivoting
function gauss_elimination(a, b, n)
    for j = 1:n-1
        partial_pivoting(a, b, n, j)
        for k = j+1:n
            m = a[k, j] / a[j, j]
            a[k, j] = m 
            a[k, j+1:n] = a[k, j+1:n] - m*a[j, j+1:n]
            b[k] = b[k] - m*b[j]
        end
    end
    return
end

# Back-Substitution
function back_substitution(a, b, n)
    for j = n:-1:1
        b[j] = (b[j] - dot(a[j, j+1:n], b[j+1:n])) / a[j, j]
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

print("Dimension: n = ")
n = parse(UInt64, readline())
n = n-1
err = 5.e-13

# Initialize A, b: Using tridiagonal matrix 
a = zeros(n, n) 
for j = 1:n
    a[j, j] = 2
end
for j = 1:n-1
    a[j, j+1] = -1
    a[j+1, j] = -1
    a[j+1, 1] = (j+1)^3
    a[1, j+1] = -(j+1)^2
    a[n, j+1] = -1/(j+1)
end
b = zeros(n)
b[1] = 1

# Copy A, b for check
savea = copy(a) 
saveb = copy(b)

gauss_elimination(a, b, n)
back_substitution(a, b, n)
check_ans(savea, saveb, b, n, err)

println("Solution: x = " )
println(b)