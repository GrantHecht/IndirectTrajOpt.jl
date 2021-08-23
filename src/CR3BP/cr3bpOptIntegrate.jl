
function cr3bpOptIntegrate(y0, tspan, ps::CR3BPIndirectParams)

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
end