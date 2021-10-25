module figures

using Plots
using PyPlot
using LaTeXStrings
using Random
using Distributions
using LinearAlgebra
using DifferentialEquations

function contour_plot(f, xrange, yrange; resolution=100)
    local n = resolution
    x = range(xrange[1],stop=xrange[2],length=n)
    y = range(yrange[1],stop=yrange[2],length=n)

    xgrid = repeat(x',n,1)
    ygrid = repeat(y,1,n)
    z = zeros(n,n)

    for i in 1:n
        for j in 1:n
            z[i:i,j:j] .= f([x[j], y[i]])
        end
    end
    p = Plots.contour(xgrid, ygrid, z, colors="black", linewidth=1.0, colorbar=false, aspect_ratio=1)
    
    return p
end


function plot_logistic_regression_data(xrange, β0, n)
    # We will need the sigmoid function
    sigmoid(x) = 1.0/(1.0+exp(-x))
  
    # Generate some data
    Random.seed!(2)
    X = zeros(n,2)
    X[:,1] .= 1.0
    X[:,2] = rand(Uniform(xrange[1], xrange[2]), n)
    Y = [rand(Bernoulli(sigmoid(β0⋅X[i,:]))) for i in 1:n]
    
    p = Plots.scatter(X[:,2], Y, label="", xlabel="X : hours", ylabel="Y : eye strain (yes/no)", yrange=(-0.1,1.25))
    
    # Overlay the sigmoid curve
    xgrid = range(xrange[1],stop=xrange[2],length=n)
    y = [sigmoid(β0[1]+β0[2]*x) for x in xgrid]
    Plots.plot!(xgrid, y, color=:green, label=L"\mathbb{P}(Y\,|\,X)", linewidth=2)
    return p
end

# Slides 7
function plot_bracket_root()
    fig, ax = subplots(1, 1, figsize=(5,5))
    ax.set_ylim([-2.0, 2.0])
    ax.set_xlim([-2.5, 2.5])
    #ax.set_xticks([])
    #ax.set_yticks([])
    x = -2:0.01:2
    f(x) = 0.5*(x+1.0)^2.0 -1.75
    y = f.(x)
    ax.plot(x, y)
    X = [-1.5, -0.5, 1.25]
    Y = [0.0, 0.0 , 0.0]
    FX = f.(X)
    col = [:blue, :green , :blue]
    ax.scatter(X,Y, s=30, marker="o", c=col)
    ax.scatter(X,FX, s=30, marker="o", c=col)
    ax.axhline(color="black")
    # Add labels
    xoffset = -0.01
    yoffset = -0.25
    ann = [
        (X[1] + xoffset, yoffset,"a"), 
        (X[2] + xoffset, yoffset, "c"), 
        (X[3] + xoffset, yoffset, "b"), 
    ]
    fs = 15
    for item in ann
        ax.annotate(item[3], xy=(item[1], item[2]), xycoords="data", textcoords="data", fontsize=fs)
    end
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    
    return fig, ax
end


function plot_Euler(sol1, sol2, λ, u0)
    t=0.0:0.01:sol1.t[end]
    y= u0*exp.(-λ *t)
    p = Plots.plot(sol1.t, sol1.u, marker=:circle, label="Forward Euler", xlabel="t", ylabel="u(t)")
    Plots.plot!(sol2.t, sol2.u, marker=:circle, label="Implicit Euler")
    Plots.plot!(t, y, label =L"u_0\,exp(-\lambda\,t)", linewidth=2)
    return p
end

function rel_osc!(du,u,p,t)
    x,y = u
    μ = p[1]
    du[1] = μ *(y - ((x^3.0)/3.0 - x ))
    du[2] = - x / μ
end

function plot_relaxation_oscillator(μval)
    u0 = [-0.1,-0.1]
    tspan = (0.0,25.0*μval)
    p = [μval]
    prob = ODEProblem(rel_osc!,u0,tspan,p)
    
    sol=solve(prob)
    
    p1 = Plots.plot(sol,vars=(1,2), label="(x(t), y(t))")
    
    p2 = Plots.plot(sol, vars=(0,1), label="x(t)")
    Plots.plot!(sol,vars=(0,2), label="y(t)")
    
    return p1, p2
end

end