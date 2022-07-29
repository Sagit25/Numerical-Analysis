# Jacobi iteration method
# Solve Ax = b

using OffsetArrays

print("Dimension: n = ")
n = parse(Int, readline())
println("Iteration: max_iter = ")
max_iter = parse(Int, readline())

# vectors to be used
x=OffsetArray(zeros(n+1),0:n)
xn=OffsetArray(zeros(n+1),0:n)

b=zeros(n-1); b[1]=n
res=zeros(n-1)

for iter = 1: max_iter
    for j = 1: n-1
    	xn[j] = (b[j] + x[j-1] + x[j+1] )/2
    end
    for j = 1: n-1 # residual xn for the error calculation
        res[j] = b[j]-(-xn[j-1]+2*xn[j]-xn[j+1]) 
    end
    global x = copy(xn) # update for the next iteration
    # do not use "x=xn", or "global x = xn"
    if iter%10 == 0
	l2norm=res'*res; l2norm = sqrt(l2norm)
	@printf("iter = %8d ;   l2norm of residual = %25.15e\n",iter, l2norm)
	if l2norm < 1.e-8
	    println( "The Jacobi iteration converged at iter < ", iter)
	    println("The solution x[0:n] is given as follows:")
	    if n==10
	          @printf("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f
		  %15.5f %15.5f %15.5f %15.5f\n", xn[0],xn[1],xn[2],xn[3],
		  xn[4],xn[5],xn[6],xn[7],xn[8],xn[9],xn[10])
	    end
#	    for j=0:10:n
#	          @printf("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f
#		  %15.5f %15.5f %15.5f %15.5f\n", xn[j:min(j+9,n)])
#            end
	    exit()
	end
    end
end
