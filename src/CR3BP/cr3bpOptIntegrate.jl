
function cr3bpOptIntegrate(y0, tspan, ps::AbstractCR3BPIndirectParams; 
    copyParams = false, termCallbacks = false, inPlace = false)

   # Instantiate problem 
   prob = createCR3BPODEProb(y0, tspan, ps; 
    copyParams = copyParams, termCallbacks = termCallbacks, inPlace = inPlace) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-12,
        abstol = 1e-12,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    # Compute time to desired final time if terminated
    timeToFinalTime = 0.0
    if termCallbacks && sol.retcode != :Success
        timeToFinalTime = tspan[2] - sol.t[end]
    end

    # Return final states and co-states
    return sol[end], timeToFinalTime
end

function cr3bpOptIntegrate(y0, tspan, ps::AbstractCR3BPIndirectParams, flag::SingleOutput; 
    copyParams = false, termCallbacks = false, inPlace = false)

   # Instantiate problem 
   prob = createCR3BPODEProb(y0, tspan, ps; 
    copyParams = copyParams, termCallbacks = termCallbacks, inPlace = inPlace) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-12,
        abstol = 1e-12,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    # Return final states and co-states
    return sol[end]
end

function cr3bpOptIntegrate(y0, tspan, ϵ, ps::AbstractCR3BPIndirectParams, flag::SingleOutput; 
    copyParams = false, termCallbacks = false, inPlace = false)

    # Instantiate problem 
    prob = createCR3BPODEProb(y0, tspan, ps; copyParams = copyParams, 
        termCallbacks = termCallbacks, inPlace = inPlace, ϵ = ϵ) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-12,
        abstol = 1e-12,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    # Return final states and co-states
    return sol[end]
end

function cr3bpOptWithSTMIntegrate(z0, tspan, ps::CR3BPIndirectWithSTMParams; copyParams = false)

   # Instantiate problem 
   prob = createCR3BPODEWithSTMProb(z0, tspan, ps; copyParams) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-12,
        abstol = 1e-12,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e10
        )

    # Return final states and co-states
    return sol[end]
end

function cr3bpOptWithSTMIntegrate(z0, tspan, ϵ, ps::CR3BPIndirectWithSTMParams; copyParams = false)

    # Instantiate problem 
    prob = createCR3BPODEWithSTMProb(z0, tspan, ps; copyParams = copyParams, ϵ = ϵ) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-12,
        abstol = 1e-12,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e10
        )

    # Return final states and co-states
    return sol[end]
end
