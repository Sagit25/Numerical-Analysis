# Gauss Jordan with partial pivoting
# Compute A^(-1)

using LinearAlgebra

# Swap columns or rows
function swap(a, j, k, op)
    if op == 1
        for n in axes(a, 1)
            a[n, j], a[n, k] = a[n, k], a[n, j]
        end
    elseif op == 2
        for n in axes(a, 2)
            a[j, n], a[k, n] = a[k, n], a[j, n]
        end
    end
    return
end

# Gauss-Jordan with partial pivoting
function gauss_jordan(a, perm, n, p)
    for j = 1:n
        if p
            js = argmax(abs.(a[j:n, j]))
            jnew = j+js-1
            perm[j] = jnew
            swap(a, j, jnew, 2)
        end
        m = 1/a[j, j]; a[j, j] = m
        col = m*a[:, j]; row = a[j, :]
        a = a-col*row'
        a[:, j] = col; a[j, :] = -m*row; a[j, j] = m
    end
    if p
        for j = n-1:-1:1
            if j != perm[j]
                swap(a, j, perm[j], 1)
            end
        end
    end
    return
end

function check_ans(a, savea, err)
    diag = zeros(n, n)
    diag = a*savea
    for j = 1:n
        for k = 1:n
            if j == k 
                if abs(diag[k, j]-1) > err
                    println("Inverse is not obtained")
                    exit()
                end
            else
                if abs(diag[k, j]) > err
                    println("Inverse is not obtained")
                    exit()
                end
            end
        end
    end
    return
end

print("Dimension: n = ")
n = parse(UInt64, readline())
print("Partial pivoting (1: use, 0: don't use): p = ")
p = parse(Bool, readline())
err = 1.e-5
    
# Initialize A, permutation matrix
a = zeros(n, n) 
for j = 1:n
    a[j, j] = 4
end
for j = 1:n-1
    a[j, j+1] = -1
    a[j+1, j] = -1
    a[min(j+3, n), j] = -5
    a[j, min(j+2, n)] = -7
end
perm = collect(1:n)
    
# Copy A for check
savea = copy(a)

gauss_jordan(a, perm, n, p)
# check_ans(a, savea, err)
    
println("Original matrix A = ")
println(savea)
println("Inverse of A = ")
println(a)
println(a*savea)
println("Condition number of A in l2-norm = ")
println(norm(a, 2)*norm(savea, p))