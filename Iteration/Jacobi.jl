# Jacobi iteration method
# Solve Ax = b (A: tridiagonal matrix)

using OffsetArrays

function jacobi_iteration(x, xn, b, max_iter)
    for iter = 1 : max_iter
        for j = 1 : n-1
            xn[j] = (b[j]+x[j-1]+x[j+1])/2
        end
        for j = 1 : n-1 
            res[j] = b[j]-(-xn[j-1]+2*xn[j]-xn[j+1]) 
        end
        global x = copy(xn) 
        if iter%10 == 0
            l2norm = sqrt(res'*res)
            println("iter = ", iter, " l2norm of residual = ", l2norm)
            if l2norm < 1.e-8
                println("Jacobi iteration converged at iter < ", iter)
                println("Solution: xn = ", xn)
                break
            end
        end
    end
end

print("Dimension: n = ")
n = parse(Int, readline())
println("Iteration: max_iter = ")
max_iter = parse(Int, readline())

x = OffsetArray(zeros(n+1),0:n)
xn = OffsetArray(zeros(n+1),0:n)

b = zeros(n-1)
b[1] = n
res = zeros(n-1)

jacobi_iteration(x, xn, b, max_iter)
