abstract type AbstractIndirectTrajOptimizer end

struct IndirectTrajOptimizer{IPT,CSIT,ST} <: AbstractIndirectOptimizationProblem
    # Indirect Trajectory Optimization problem 
    prob::IPT

    # Indirect Co-State Initializer 
    csInit::CSIT

    # Indirect solver 
    solver::ST

    # Data Output Manager
    writingData::Bool
    dataOutputManager::DataOutputManager

    # Meta data (May want to add more here) 
    solMethod::Symbol
    initCost::Symbol 
    initOptimizer::Symbol
    weights::Vector{Float64}
    minNeighborhoodFrac::Float64
    UBs::Vector{Float64}
    LBs::Vector{Float64}
    iUBs::Vector{Float64}
    iLBs::Vector{Float64}
    maxIters::Int 
    funcTol::Float64
    maxStallIters::Float64
    maxStallTime::Float64
    maxTime::Float64
    useParallel::Bool 
end

function IndirectTrajOptimizer(prob; solutionMethod = :FSS, nSeg = 4, initCostFunc = :WSS, 
    initOptimizer = :PSO, numParticles = 500, numSwarms = 4, swarmInitMethod = :Uniform, 
    UBs = [100, 100, 100, 50, 50, 50, 50], LBs = [-100, -100, -100, -50, -50, -50, -50],
    iUBs = nothing, iLBs = nothing, weights = [10, 10, 10, 1, 1, 1, 1], MFD = 2.8e-6, display = true,
    displayInterval = 1, maxIters = 1000, funcTol = 1e-6, maxStallIters = 25, maxStallTime = 500,
    maxTime = 1800, useParallel = true, homotopyParamVec = nothing, dataFolder = nothing,
    minNeighborhoodFraction = 0.25,writeData = false, initCallback = nothing)

    # Check that homotopy parameter vector has been set if homotopy is used 
    if prob.homotopy && homotopyParamVec === nothing
        throw(ArgumentError("Problem specified to use homotopy continuation but continuation parameters not provided."))
    end

    # Initialize initializer and solver 
    if solutionMethod == :FSS
        # Initialize co-state initializer
        csInit = FSSCoStateInitializer(prob.initBVPFunc, prob.tspan, prob.iConds, prob.fConds;
            costFunc = initCostFunc, optimizer = initOptimizer, numParticles = numParticles,
            numSwarms = numSwarms, initMethod = swarmInitMethod, UBs = UBs, LBs = LBs,
            iUBs = iUBs, iLBs = iLBs, weights = weights, display = display,
            displayInterval = displayInterval, maxIters = maxIters, funcTol = funcTol, 
            MFD = MFD, minNeighborhoodFraction = minNeighborhoodFraction,
            maxStallIters = maxStallIters, maxStallTime = maxStallTime, maxTime = maxTime,
            useParallel = useParallel)

        # Initialize solver
        solver = FSSSolver(zeros(7), tspan, prob.iConds, prob.fConds, prob.BVPFunc, prob.BVPWithSTMFunc;
            homotopy = prob.homotopy, homotopyParamVec = homotopyParamVec)
    elseif solutionMethod == :FMS
        # Initialize co-state initializer
        csInit = FSSCoStateInitializer(prob.initBVPFunc, prob.tspan, prob.iConds, prob.fConds;
            costFunc = initCostFunc, optimizer = initOptimizer, numParticles = numParticles,
            numSwarms = numSwarms, initMethod = swarmInitMethod, UBs = UBs, LBs = LBs,
            iUBs = iUBs, iLBs = iLBs, weights = weights, display = display,
            displayInterval = displayInterval, maxIters = maxIters, funcTol = funcTol, 
            MFD = MFD, minNeighborhoodFraction = minNeighborhoodFraction,
            maxStallIters = maxStallIters, maxStallTime = maxStallTime, maxTime = maxTime,
            useParallel = useParallel)

        # Initialize solver
        solver = FMSSolver(tspan, prob.iConds, prob.fConds, prob.BVPFunc, prob.BVPWithSTMFunc;
            homotopy = prob.homotopy, homotopyParamVec = homotopyParamVec, nSeg = nSeg)
    else
        throw(ArgumentError("Only forward single and multiple shooting are implemented."))
    end

    # Initialize DataOutputManager 
    if dataFolder === nothing 
        dataFolder = joinpath(pwd(),"data")
    end
    dataOutputManager = DataOutputManager(dataFolder)

    if iUBs === nothing 
        iUBs = Vector{Float64}(undef, 0)
    end
    if iLBs === nothing 
        iLBs = Vector{Float64}(undef, 0)
    end

    IndirectTrajOptimizer{typeof(prob), typeof(csInit), typeof(solver)}(prob, csInit, solver, 
        writeData, dataOutputManager, solutionMethod, initCostFunc, initOptimizer, weights,
        minNeighborhoodFraction, UBs, LBs, iUBs, iLBs, maxIters, funcTol, maxStallIters, 
        maxStallTime, maxTime, useParallel)
end

function initialize!(ito::IndirectTrajOptimizer)
    # Initialize co-states
    initialize!(ito.csInit)

    # Re-initialize solvers data manager
    initializeData!(ito.solver, GetInitializedCostates(ito.csInit))

    # Write data if desired 
    if ito.writingData
        writeData(ito.dataOutputManager, ito)
    end

    return nothing
end

function solve!(ito::IndirectTrajOptimizer; factor = 3.0, ftol = 1e-8, showTrace = true, convergenceAttempts = 4)

    # Solve 
    solve!(ito.solver; factor = factor, ftol = ftol, showTrace = showTrace, convergenceAttempts = convergenceAttempts)

    # Write data if desired 
    if ito.writingData
        writeData(ito.dataOutputManager, ito)
    end

    return nothing
end

function tSolve!(itoVec; factor = 3.0, ftol = 1e-8, showTrace = true, convergenceAttempts = 4)
    p = Progress(length(itoVec), 1, "Solving BVPs: ")
    Threads.@threads for i in eachindex(itoVec)
        solve!(itoVec[i]; factor = factor, ftol = ftol, showTrace = showTrace, convergenceAttempts = convergenceAttempts)
        next!(p)
    end
    return nothing
end
