# Function writes all data from the IndirectTrajectoryOptimizer to 
# a .txt file and an HDL5 binary file (using JLD2.jl).
#
# Called after each action is performed by IndirectTrajOptimizer to 
# "hopefully" ensure no data loss.
function writeData(dom::DataOutputManager, ito::AbstractIndirectTrajOptimizer)

    # Initialize file name if necessary. If file name is not initialized
    # directories may not be either so initialize directories if necessary.
    if dom.fnameInitialized == false
        # Check that directories exists and create it if not 
        if !isdir(dom.baseFolder)
            mkpath(dom.baseFolder)
        end
        if !isdir(joinpath(dom.baseFolder, "textData"))
            mkdir(joinpath(dom.baseFolder, "textData"))
        end
        if !isdir(joinpath(dom.baseFolder, "binaryData"))
            mkdir(joinpath(dom.baseFolder, "binaryData"))
        end

        # Initialize fname to ensure no data is overwritten
        dom.fname = "results_#1"
        if isfile(joinpath(dom.baseFolder, "textData", dom.fname*".txt"))
            safe = false
            cnt  = 1
            while !safe 
                cnt += 1
                dom.fname = "results_#" * string(cnt)
                if !isfile(joinpath(dom.baseFolder, "textData", dom.fname*".txt"))
                    # Create new text file 
                    touch(joinpath(dom.baseFolder, "textData", dom.fname*".txt"))
 
                    # Set flags
                    safe = true
                    dom.fnameInitialized = true
               end
            end
        end
    end

    # Write data
    writeBinaryData(dom, ito)
    writeTextData(dom, ito)

    return nothing
end

function writeBinaryData(dom::DataOutputManager, ito::AbstractIndirectTrajOptimizer)
    bson(joinpath(dom.baseFolder, "binaryData", dom.fname*".bson"), Dict(:data=>ito))
    return nothing
end

function writeTextData(dom::DataOutputManager, ito::IndirectTrajOptimizer)
    # Open file
    fid = open(joinpath(dom.baseFolder, "textData", dom.fname*".txt"), "w")

    # Write meta data
    println(fid, "# META DATA")
    println(fid, "File Format Version:\t" * dom.textFormatVersion)
    println(fid, "Solution Method:\t\t" * string(ito.solMethod))
    println(fid, "Initialization Cost:\t" * string(ito.initCost))
    println(fid, "Heuristic Optimizer:\t" * string(ito.initOptimizer))
    if ito.initOptimizer == :PSO
        numParticles = length(ito.csInit.ho.swarm)
        println(fid, "Number of Particles:\t" * string(numParticles))
    elseif ito.initOptimizer == :MS_PSO
        numParticles = length(ito.csInit.ho.swarmVec[1])
        numSwarms    = length(ito.csInit.ho.swarmVec)
        println(fid, "Number of Particles:\t" * string(numParticles))
        println(fid, "Number of Swarms:\t\t" * string(numSwarms))
    end                                                           
    println(fid, "Using Homotopy:\t\t\t" * (ito.prob.homotopy ? "Yes" : "No"))
    println(fid, "Min. Neighborhood Frac:\t" * string(ito.minNeighborhoodFrac))
    if ito.initCost == :WSS || ito.initCost == :WSSWM
        lineStr = "Cost Funciton Weights:\t["
        for i in 1:length(ito.weights)
            if i == length(ito.weights)
                lineStr *= string(ito.weights[i])
            else
                lineStr *= string(ito.weights[i]) * ", "
            end
        end
        lineStr *= "]"
        println(fid, lineStr)
    end
    lineStr = "Initialization UBs:\t\t["
    for i in 1:length(ito.UBs)
        if i == length(ito.UBs)
            lineStr *= string(ito.UBs[i])
        else
            lineStr *= string(ito.UBs[i]) * ", "
        end
    end
    lineStr *= "]"
    println(fid, lineStr)
    lineStr = "Initialization LBs:\t\t["
    for i in 1:length(ito.LBs)
        if i == length(ito.LBs)
            lineStr *= string(ito.LBs[i])
        else
            lineStr *= string(ito.LBs[i]) * ", "
        end
    end
    println(fid, lineStr)
    if length(ito.iUBs) != 0
        lineStr = "Initialization iUBs:\t["
        for i in 1:length(ito.iUBs)
            if i == length(ito.iUBs)
                lineStr *= string(ito.iUBs[i])
            else
                lineStr *= string(ito.iUBs[i]) * ", "
            end
        end
        lineStr *= "]"
        println(fid, lineStr)
    end
    if length(ito.iLBs) != 0
        lineStr = "Initialization iLBs:\t["
        for i in 1:length(ito.iLBs)
            if i == length(ito.iLBs)
                lineStr *= string(ito.iLBs[i])
            else
                lineStr *= string(ito.iLBs[i]) * ", "
            end
        end
        lineStr *= "]"
        println(fid, lineStr)
    end
    println(fid, "Initialization Func. Tol.:\t\t\t\t"*string(ito.funcTol))
    println(fid, "Initialization Max Iterations:\t\t\t"*string(ito.maxIters))
    println(fid, "Initialization Max Stall Iterations:\t"*string(ito.maxStallIters))
    println(fid, "Initialization Max Time:\t\t\t\t"*string(ito.maxTime)*" sec")
    println(fid, "Initialization Max Stall Time:\t\t\t"*string(ito.maxStallTime)*" sec")
    println(fid, "Used Parallel:\t\t\t\t\t\t\t" * (ito.useParallel ? "Yes" : "No"))

    println(fid, "# END META DATA"); println(fid, "")

    # Write convergence data
    println(fid, "# CONVERGENCE DATA")
    if !(ito.csInit.ho isa Symbol)
        println(fid, "Time to Initialize:\t\t\t\t" * string(ito.csInit.ho.results.time) * " sec")
        println(fid, "Initialization Iterations:\t\t" * string(ito.csInit.ho.results.iters))
        if ito.initOptimizer == :PSO || ito.initOptimizer == :MS_PSO
            println(fid, "Initialization Func. Evals.:\t" * 
                string(numParticles*ito.csInit.ho.results.iters))
        end
        println(fid, "Heuristic Opj. Function:\t\t" * string(ito.csInit.ho.results.fbest))
    end
    println(fid, "Initial Guess Converged:\t\t" * 
        (GetInitialGuessConverged(ito.solver) ? "Yes" : "No"))
    if ito.prob.homotopy 
        println(fid, "Homotopy Converged:\t\t\t\t" * 
        (GetHomotopyConverged(ito.solver) ? "Yes" : "No"))
    end
    println(fid, "# END CONVERGENCE DATA"); println(fid, "")

    # Write co-state data
    println(fid, "# COSTATE DATA")
    println(fid, "Initialized Co-States:")
    initCSVec = GetInitializedCostates(ito.csInit)
    linestr = ""
    for i in 1:length(initCSVec)
        linestr *= string(initCSVec[i]) * "\t"
    end
    println(fid, linestr)
    if ito.prob.homotopy 
        solVec = GetHomotopySolutionVector(ito.solver)
        ϵs = GetHomotopyParams(ito.solver)
        cflags = GetHomotopyConvergenceFlags(ito.solver)
        for i in 1:length(solVec)
            println(fid, "param: " * string(ϵs[i]) * " converged: " * (cflags[i] ? "Yes" : "No"))
            linestr = ""
            for j in 1:length(solVec[i])
                linestr *= string(solVec[i][j]) * "\t" 
            end
            println(fid, linestr)
        end
    else
        println(fid, "converged " * (GetInitialGuessConverged(ito.solver) ? "Yes" : "No"))
        linestr = ""
        sol = GetSolution(ito.solver)
        for i in 1:length(sol)
            linestr *= string(sol[i]) * "\t"
        end
        println(fid, linestr)
    end
    println(fid, "# END COSTATE DATA")

    # Close file
    close(fid)
    return nothing
end

function writeTextData(dom::DataOutputManager, ito::InitializedIndirectTrajOptimizer)
    # Open file
    fid = open(joinpath(dom.baseFolder, "textData", dom.fname*".txt"), "w")

    # Write meta data
    println(fid, "# META DATA")
    println(fid, "File Format Version:\t" * dom.textFormatVersion * "i")
    println(fid, "Solution Method:\t\t" * string(ito.solMethod))
    println(fid, "Using Homotopy:\t\t\t" * (ito.prob.homotopy ? "Yes" : "No"))
    println(fid, "# END META DATA"); println(fid, "")

    # Write convergence data
    println(fid, "# CONVERGENCE DATA")
    if ito.time > 0.0
        println(fid, "Time to Initialize:\t\t\t\t" * string(ito.time) * " sec")
    end
    if ito.iters > 0
        println(fid, "Initialization Iterations:\t\t" * string(ito.iters))
    end
    if ito.fevals > 0
        println(fid, "Initialization Func. Evals.:\t" * 
            string(ito.fevals))
    end
    if ito.fval >= 0.0
        println(fid, "Heuristic Opj. Function:\t\t" * string(ito.fval))
    end
    println(fid, "Initial Guess Converged:\t\t" * 
        (GetInitialGuessConverged(ito.solver) ? "Yes" : "No"))
    if ito.prob.homotopy 
        println(fid, "Homotopy Converged:\t\t\t\t" * 
        (GetHomotopyConverged(ito.solver) ? "Yes" : "No"))
    end
    println(fid, "# END CONVERGENCE DATA"); println(fid, "")

    # Write co-state data
    println(fid, "# COSTATE DATA")
    println(fid, "Initialized Co-States:")
    initCSVec = ito.λi
    linestr = ""
    for i in 1:length(initCSVec)
        linestr *= string(initCSVec[i]) * "\t"
    end
    println(fid, linestr)
    if ito.prob.homotopy 
        solVec = GetHomotopySolutionVector(ito.solver)
        ϵs = GetHomotopyParams(ito.solver)
        cflags = GetHomotopyConvergenceFlags(ito.solver)
        for i in 1:length(solVec)
            println(fid, "param: " * string(ϵs[i]) * " converged: " * (cflags[i] ? "Yes" : "No"))
            linestr = ""
            for j in 1:length(solVec[i])
                linestr *= string(solVec[i][j]) * "\t" 
            end
            println(fid, linestr)
        end
    else
        println(fid, "converged " * (GetInitialGuessConverged(ito.solver) ? "Yes" : "No"))
        linestr = ""
        sol = GetSolution(ito.solver)
        for i in 1:length(sol)
            linestr *= string(sol[i]) * "\t"
        end
        println(fid, linestr)
    end
    println(fid, "# END COSTATE DATA")

    # Close file
    close(fid)
    return nothing
end
