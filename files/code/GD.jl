module GD

using LinearAlgebra

include("gss.jl")

# Gradient descent with full line minimisation


# Gradient descent with a pre-defined learning rate, η(k)
function gradient_descent(f, df, x0, ϵtol, maxiter, η)
    g = df(x0)
    n = 0 
    x = x0
    dx = - g / norm(g)

    steps = Int[]
    xposition = Float64[]
    yposition = Float64[]

    while n<maxiter && norm(g)>ϵtol
        append!(steps,n); append!(xposition, x[1]); append!(yposition, x[2]);
        # Step size is externally prescribed
        s = η(n+1)
        x = x + s*dx
        g = df(x)
        dx = -g/norm(g)
        n+=1
    end
    return xposition, yposition
end

# Vanilla stochastic gradient descent
function stochastic_gradient_descent(f, df, x0, maxiter, η, datasize, batchsize)
    g = df(x0,1)
    n = 0 
    x = x0
    dx = - g

    steps = Int[]
    xposition = Float64[]
    yposition = Float64[]

    while n<maxiter
        append!(steps,n); append!(xposition, x[1]); append!(yposition, x[2]);
        # Step size is externally prescribed
        s = η(n+1)
        x = x + s*dx
        i = 1+(n%(datasize-1))
        g = df(x, i)
        dx = -g
        n+=1
    end
    return xposition, yposition
end

end
