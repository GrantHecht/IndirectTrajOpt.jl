
struct IndirectTrajOptimizer{IPT,CSIT,ST}
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
end

function IndirectTrajOptimizer(prob; solutionMethod = :FSS, initCostFunc = :WSS, 
    initOptimizer = :PSO, numParticles = 500, numSwarms = 4, swarmInitMethod = :Uniform, 
    UBs = [100, 100, 100, 50, 50, 50, 50], LBs = [-100, -100, -100, -50, -50, -50, -50],
    iUBs = nothing, iLBs = nothing, weights = [10, 10, 10, 1, 1, 1, 1], display = true,
    displayInterval = 1, maxIters = 1000, funcTol = 1e-6, maxStallIters = 25, maxStallTime = 500,
    maxTime = 1800, useParallel = true, homotopyParamVec = nothing, dataFolder = nothing,
    writeData = false)

    # Check that homotopy parameter vector has been set if homotopy is used 
    if prob.homotopy && homotopyParamVec === nothing
        throw(ArgumentError("Problem specified to use homotopy continuation but continuation parameters not provided."))
    end

    # Initialize initializer and solver 
    if solutionMethod == :FSS

        # Initialize co-state initializer
        csInit = FSSCoStateInitializer(prob.initBVPFunc, prob.iConds, prob.fConds;
            costFunc = initCostFunc, optimizer = initOptimizer, numParticles = numParticles,
            numSwarms = numSwarms, initMethod = swarmInitMethod, UBs = UBs, LBs = LBs,
            iUBs = iUBs, iLBs = iLBs, weights = weights, display = display,
            displayInterval = displayInterval, maxIters = maxIters, funcTol, 
            maxStallIters = maxStallIters, maxStallTime = maxStallTime, maxTime = maxTime,
            useParallel = useParallel)

        # Initialize solver
        solver = FSSSolver(zeros(7), prob.iConds, prob.fConds, prob.BVPFunc, prob.BVPWithSTMFunc;
            homotopy = prob.homotopy, homotopyParamVec = homotopyParamVec)
    else
        throw(ArgumentError("Only forward sigle shooting is implemented."))
    end

    # Initialize DataOutputManager 
    if dataFolder === nothing 
        dataFolder = joinpath(pwd(),"data")
    end
    dataOutputManager = DataOutputManager(dataFolder)

    IndirectTrajOptimizer{typeof(prob), typeof(csInit), typeof(solver)}(prob, csInit, solver, 
        writeData, dataOutputManager, solutionMethod, initCostFunc, initOptimizer, weights)
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

function solve!(ito::IndirectTrajOptimizer)
    solve!(ito.solver)

    # Write data if desired 
    if ito.writingData
        writeData(ito.dataOutputManager, ito)
    end

    return nothing
end