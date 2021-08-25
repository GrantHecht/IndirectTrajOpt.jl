
function cr3bpOptIntegrate(y0, tspan, ps::AbstractCR3BPIndirectParams)

   # Instantiate problem 
   prob = createCR3BPODEProb(y0, tspan, ps) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    # Return final states and co-states
    return sol[end]
end

function cr3bpOptWithSTMIntegrate(z0, tspan, ps::CR3BPIndirectWithSTMParams)

   # Instantiate problem 
   prob = createCR3BPODEWithSTMProb(z0, tspan, ps) 

    # Solve ode 
    sol = solve(
        prob,
        TsitPap8(),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    # Return final states and co-states
    return sol[end]
end