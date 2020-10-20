module gss

# Function to implement golden section search
function minimise(f, triple, tol; verbose=false)
    if verbose
        # Some arrays to store intermediate values for plotting
        steps = Int64[]
        position = Float64[]
        width = Float64[]
    end
    
    # Initialise
    # It is a good idea to ensure that the initial bracketing triple has the right
    # ratio in case the function caller forgets to do this.
    ϕ = 0.5*(sqrt(5.0)-1.0)
    a = triple[1]
    b = triple[2]
    c = b - (b-a)*ϕ
    (fa, fc, fb) = (f(a), f(c), f(b))
     # Check that this is indeed a bracketing triple
    if !(fa > fc && fb > fc && b > a)
        error(triple, " does not bracket a minimum: values are ", (fa, fb, fc))
    end
        
    n = 0 
    
    while(b-a > tol && n <30)
        x = a + b -c
        fx = f(x)
        
        if verbose
            append!(steps,n)
            append!(position, x)
            append!(width, b-a)
        end
        
        if (c-a < b-c)
            x1 = c; f1 = fc
            x2 = x; f2 = fx
        else
            x1=x; f1 = fx
            x2=c; f2 = fc
        end
        
        if (f1 < f2)
            c = x1
            b = x2
            fc = f1
            fb = f2
        else
            a = x1
            c = x2
            fa = f1
            fc = f2
        end
        n+=1
        
    end
    
    if verbose
        return steps, position, width
    else
        return c
    end
end

end