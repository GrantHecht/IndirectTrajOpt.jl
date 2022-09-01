
function integrate(y0, tspan, ps::AbstractCR3BPIndirectParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::InitializationFlag, homotopyFlag::HomotopyFlag; 
    copyParams = false, termCallbacks = false, inPlace = false)

   # Instantiate problem 
   prob = createCR3BPODEProb(y0, tspan, ps, integrationTypeFlag, homotopyFlag; 
        copyParams = copyParams, termCallbacks = termCallbacks, inPlace = inPlace) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e7
        )

    # Compute time to desired final time if terminated
    timeToFinalTime = 0.0
    if termCallbacks && sol.retcode != :Success
        timeToFinalTime = tspan[2] - sol.t[end]
    end

    # Return final states and co-states
    return sol[end], timeToFinalTime
end

function integrate(y0, tspan, ps::CR3BPIndirectWithSTMParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::InitializationWithJacobianRankPenalty, homotopyFlag::HomotopyFlag; 
    copyParams = false, termCallbacks = false)

    # Create initial state with STM
    z0   = zeros(210)
    @inbounds for i in 1:14
        # State/Co-State Entries
        z0[i] = y0[i]

        # STM Entries
        idx     = 14 + (i - 1)*14 + i
        z0[idx] = 1.0
    end

    # Instantiate problem 
    prob = createCR3BPODEWithSTMProb(z0, tspan, ps, homotopyFlag; copyParams = copyParams) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy=false),
        reltol = 1e-10,
        abstol = 1e-10,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6,
        verbose = false
    )

    # Process output solution to get final state and Jacobian rank penalty
    yf  = sol[end][1:14]

    # Not efficient but doing this for now!
    stm = reshape(sol[end][15:end], (14, 14))
    jac = stm[[1,2,3,4,5,6,14], 8:14]

    # Compute time to desired final time if terminated
    timeToFinalTime = 0.0
    if sol.retcode != :Success
        timeToFinalTime = tspan[2] - sol.t[end]
    end
    
    # Compute penalty
    penalty = timeToFinalTime
    if rank(jac) < 7
        penalty += 1.0
    end

    # Return final states and co-states
    return yf, penalty
end

function integrate(y0, tspan, ps::AbstractCR3BPIndirectParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::Solving, homotopyFlag::HomotopyFlag; 
    copyParams = false, termCallbacks = false, inPlace = false)

    # Instantiate problem 
    prob = createCR3BPODEProb(y0, tspan, ps, integrationTypeFlag, homotopyFlag; 
        copyParams = copyParams, termCallbacks = termCallbacks, inPlace = inPlace) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e7
        )

    # Return final states and co-states
    return sol[end]
end

function integrate(y0, tspan, ps::AbstractCR3BPIndirectParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::FullSolutionHistory, homotopyFlag::HomotopyFlag; 
    copyParams = false, termCallbacks = false, inPlace = false)

    # Instantiate problem 
    prob = createCR3BPODEProb(y0, tspan, ps, integrationTypeFlag, homotopyFlag; 
        copyParams = copyParams, 
    termCallbacks = termCallbacks, inPlace = inPlace, save_positions = (true, true)) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        maxiters = 1e7
        )

    # Return final states and co-states
    return sol
end

function integrateWithHomotopy(y0, tspan, ϵ, ps::AbstractCR3BPIndirectParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::Solving, homotopyFlag::HomotopyFlag; 
    copyParams = false, termCallbacks = false, inPlace = false)

    # Instantiate problem 
    prob = createCR3BPODEProb(y0, tspan, ps, integrationTypeFlag, homotopyFlag; copyParams = copyParams, 
        termCallbacks = termCallbacks, inPlace = inPlace, ϵ = ϵ) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e7
        )

    # Return final states and co-states
    return sol[end]
end

function integrate(z0, tspan, ps::CR3BPIndirectWithSTMParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::SolvingWithSTM, homotopyFlag::HomotopyFlag; 
    copyParams = false)

   # Instantiate problem 
   prob = createCR3BPODEWithSTMProb(z0, tspan, ps, homotopyFlag; copyParams) 

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e10
        )

    # Return final states and co-states
    return sol[end]
end

function integrateWithHomotopy(z0, tspan, ϵ, ps::CR3BPIndirectWithSTMParams, dynamicsFlag::CR3BP, 
    integrationTypeFlag::SolvingWithSTM, homotopyFlag::HomotopyFlag; copyParams = false)

    # Instantiate problem 
    prob = createCR3BPODEWithSTMProb(z0, tspan, ps, homotopyFlag; copyParams = copyParams, ϵ = ϵ) 
    
    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e10
        )

    # Return final states and co-states
    return sol[end]
end

function integrate(x0, tspan, scenario::String, dynamicsFlag::CR3BP, 
    integrationTypeFlag::FullSolutionHistoryNoControl; inPlace = false)

   # Instantiate problem 
   prob = createCR3BPODENoControlProb(x0, tspan, scenario; inPlace = inPlace)

    # Solve ode 
    sol = solve(
        prob,
        Vern9(lazy = false),
        reltol = 1e-14,
        abstol = 1e-14,
        maxiters = 1e7
        )

    # Return final states and co-states
    return sol
end
