
mutable struct DataOutputManager

    # Base folder location. Should be full file path and not local
    baseFolder::String

    # File name without extension
    fname::String

    # File name initialized flag
    fnameInitialized::Bool

    # Text file format version
    textFormatVersion::String

    # Constructor
    function DataOutputManager(baseFolder::String)
        new(baseFolder, "", false, "v0.1")
    end
end

# Function writes all data from the IndirectTrajectoryOptimizer to 
# a .txt file and an HDL5 binary file (using JLD2.jl).
#
# Called after each action is performed by IndirectTrajOptimizer to 
# "hopefully" ensure no data loss.
function writeData(dom::DataOutputManager, ito)

    # Initialize file name if necessary. If file name is not initialized
    # directories may not be either so initialize directories if necessary.
    if !dom.fnameInitialized
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
                    safe = true
                    dom.fnameInitialized = true
                end
            end
        end

        # Create new text file 
        touch(joinpath(dom.baseFolder, "textData", dom.fname*".txt"))
    end

    # Write data
    writeBinaryData(dom, ito)
    writeTextData(dom, ito)

    return nothing
end

function writeBinaryData(dom::DataOutputManager, ito)
    jldsave(joinpath(dom.baseFolder, "binaryData", dom.fname*".jld2"), data = ito)
    return nothing
end

function writeTextData(dom::DataOutputManager, ito)
    # Open file
    fid = open(joinpath(dom.baseFolder, "textData", dom.fname*".txt"), "w")

    # Write meta data
    println(fid, "# META DATA")
    println(fid, "File Format Version:\t" * dom.textFormatVersion)
    println(fid, "Solution Method:\t\t" * string(ito.solMethod))
    println(fid, "Initialization Cost:\t" * string(ito.initCost))
    if ito.initCost == :WSS || ito.initCost == :WSSWM
        lineStr = "Cost Funciton Weights:\t["
        for i in 1:length(ito.weights)
            lineStr *= string(ito.weights[i]) * ", "
        end
        lineStr *= "]"
        println(fid, lineStr)
    end
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
    println(fid, "# END META DATA"); println(fid, "")

    # Write convergence data
    println(fid, "# CONVERGENCE DATA")
    println(fid, "Time to Initialize:\t\t\t\t" * string(ito.csInit.ho.results.time) * " sec")
    println(fid, "Initialization Iterations:\t\t" * string(ito.csInit.ho.results.iters))
    if ito.initOptimizer == :PSO || ito.initOptimizer == :MS_PSO
        println(fid, "Initialization Func. Evals.:\t" * 
            string(numParticles*ito.csInit.ho.results.iters))
    end
    println(fid, "Heuristic Opj. Function:\t\t" * string(ito.csInit.ho.results.fbest))
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
            if cflags[i]
                for j in 1:length(solVec[i])
                    linestr *= string(solVec[i][j]) * "\t" 
                end
            end
            println(fid, linestr)
        end
    else
        println(fid, "converged " * (GetInitialGuessConverged(ito.solver) ? "Yes" : "No"))
        linestr = ""
        if GetInitialGuessConverged(ito.solver)
            sol = GetSolution(ito.solver)
            for i in 1:length(sol)
                linestr *= string(sol[i]) * "\t"
            end
        end
        println(fid, linestr)
    end
    println(fid, "# END COSTATE DATA")

    # Close file
    close(fid)
    return nothing
end