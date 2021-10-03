
# Function for reading binary data files (.jld2) into memory. Returns vector
# of IndirectTrajOptimizer(s).
#
# folder:   String containing the full file path to folder containing binary
#           data files.
#
function readBinaryData(folder::String)
    dataVec = try
        readBinaryDataConcrete(folder) 
    catch
        readBinaryDataAbstract(folder)
    end
    return dataVec
end

function readBinaryDataAbstract(folder::String)

    # Get all files in folder
    files = readdir(joinpath(folder, "binaryData"))

    # Read all files in push data to data vector
    local dataVec
    @showprogress "Reading Binary Data Files: " for i in 1:length(files)
        # Open file
        #f = jldopen(joinpath(folder, "binaryData", files[i]), "r")
        data = get(BSON.load(joinpath(folder, "binaryData", files[i]), @__MODULE__), :data, nothing)
        if data === nothing
            throw(ErrorException("Binary file was not read succesfully."))
        end

        # Grab data and push to vector
        @suppress begin
            if i == 1
                dataVec = Vector{AbstractIndirectOptimizationProblem}([data])
            else
                push!(dataVec, data)
            end
        end
    end

    return dataVec
end

function readBinaryDataConcrete(folder::String)

    # Get all files in folder
    files = readdir(joinpath(folder, "binaryData"))

    # Read all files in push data to data vector
    local dataVec
    @showprogress "Reading Binary Data Files: " for i in 1:length(files)
        # Open file
        #f = jldopen(joinpath(folder, "binaryData", files[i]), "r")
        data = get(BSON.load(joinpath(folder, "binaryData", files[i]), @__MODULE__), :data, nothing)
        if data === nothing
            throw(ErrorException("Binary file was not read succesfully."))
        end

        # Grab data and push to vector
        @suppress begin
            if i == 1
                dataVec = Vector{AbstractIndirectOptimizationProblem}([data])
            else
                push!(dataVec, data)
            end
        end
    end

    return dataVec
end

function readTextData(folder::String)
    # Get all files in folder 
    relFiles = readdir(joinpath(folder, "textData"))
    files = map(rf->joinpath(folder, "textData", rf), relFiles)

    # Get text file version number 
    ver = getTextFileVersion(files)

    # Call function to read text file with ver version number
    if ver == v"0.1"
        df = readTextDataV0_1(files)
    else
        throw(ErrorException("Reader for text file version " * string(ver) *
            " is not implemented!"))
    end
    
    return df
end 

function readTextDataV0_1(files)

    # Instantiate DataFrame
    df = DataFrame(SolutionMethod = Symbol[], InitializationCost = Symbol[],
        CostFunctionWeights = Vector{Vector{Float64}}(undef, 0), HeuristicOptimizer = Symbol[],
        NumberOfParticles = Int[], UsingHomotopy = Bool[], TimeToInitialize = Float64[],
        InitializationIterations = Int[], InitializationFunctionEvals = Int[],
        HeuristicObjectiveFunction = Float64[], InitialGuessConverged = Bool[],
        HomotopyConverged = Bool[], InitializedCoStates = Vector{Vector{Float64}}(undef, 0),
        HomotopySolutions = Vector{Vector{Vector{Float64}}}(undef, 0),
        HomotopyParams = Vector{Vector{Float64}}(undef, 0), 
        HomotopyStepConverged = Vector{Vector{Bool}}(undef, 0),
        InitializationUBs = Vector{Vector{Float64}}(undef, 0),
        InitializationLBs = Vector{Vector{Float64}}(undef, 0),
        InitializationiUBs = Vector{Vector{Float64}}(undef, 0),
        InitializationiLBs = Vector{Vector{Float64}}(undef, 0),
        MinNeighborhoodFraction = Vector{Float64}(undef, 0),
        InitializationFuncTol = Vector{Float64}(undef, 0),
        InitializationMaxIters = Vector{Int}(undef, 0),
        InitializationMaxStallIters = Vector{Int}(undef, 0),
        InitializationMaxTime = Vector{Float64}(undef, 0),
        InitializationMaxStallTime = Vector{Float64}(undef, 0))

    @showprogress "Reading Text Data Files: " for i in 1:length(files)
        # Open file
        f = open(files[i], "r")

        # Initialize local data 
        solMethod   = :NA
        initCost    = :NA 
        cFW         = Vector{Float64}(undef, 0)
        ho          = :NA
        numP        = 0
        usingH      = false
        toi         = 0.0
        initIters   = 0
        initFEvals  = 0
        hoFVal      = 0.0
        initConv    = false
        hConv       = false
        initCS      = Vector{Float64}(undef, 0)
        hSols       = Vector{Vector{Float64}}(undef, 0)
        hParams     = Vector{Float64}(undef, 0)
        hStepConv   = Vector{Bool}(undef, 0)
        ubs         = Vector{Float64}(undef, 0)
        lbs         = Vector{Float64}(undef, 0)
        iubs        = Vector{Float64}(undef, 0)
        ilbs        = Vector{Float64}(undef, 0)
        mnf         = -1.0
        funcTol     = 0.0
        maxIters    = 0
        maxSIters   = 0
        maxTime     = 0.0
        maxSTime    = 0.0

        # File location flags and line counter
        inMetaData      = false
        inConvData      = false
        inCoStateData   = false
        csDataCnt       = 1

        # Loop through lines
        for line in readlines(f)

            # Read data depending on file location
            if inMetaData
                if occursin("Solution Method:", line)
                    solMethod = Symbol(split(line, "\t"; keepempty = false)[2])
                elseif occursin("Initialization Cost:", line)
                    initCost = Symbol(split(line, "\t"; keepempty = false)[2])
                elseif occursin("Heuristic Optimizer", line)
                    ho = Symbol(split(line, "\t"; keepempty = false)[2])
                elseif occursin("Number of Particles:", line)
                    numP = parse(Int64, split(line, "\t"; keepempty = false)[2])
                elseif occursin("Using  Homotopy:", line)
                    if occursin("Yes", line)
                        usingH = true
                    end
                elseif occursin("Cost Funciton Weights:", line)
                    # Process string
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "[" => "")
                    str = replace(str, "]" => "")
                    str = replace(str, " " => "")[1:end-1]
                    strVec = split(str, ","; keepempty = false)

                    # Fill vector
                    resize!(cFW, length(strVec))
                    cFW .= map(s->parse(Float64, s), strVec)
                elseif occursin("Initialization UBs:", line)
                    # Process string
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "[" => "")
                    str = replace(str, "]" => "")
                    str = replace(str, " " => "")[1:end-1]
                    strVec = split(str, ","; keepempty = false)

                    # Fill vector
                    resize!(ubs, length(strVec))
                    ubs .= map(s->parse(Float64, s), strVec)
                elseif occursin("Initialization LBs:", line)
                    # Process string
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "[" => "")
                    str = replace(str, "]" => "")
                    str = replace(str, " " => "")[1:end-1]
                    strVec = split(str, ","; keepempty = false)

                    # Fill vector
                    resize!(lbs, length(strVec))
                    lbs .= map(s->parse(Float64, s), strVec)
                elseif occursin("Initialization iUBs:", line)
                    # Process string
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "[" => "")
                    str = replace(str, "]" => "")
                    str = replace(str, " " => "")[1:end-1]
                    strVec = split(str, ","; keepempty = false)

                    # Fill vector
                    resize!(iubs, length(strVec))
                    iubs .= map(s->parse(Float64, s), strVec)
                elseif occursin("Initialization iLBs:", line)
                    # Process string
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "[" => "")
                    str = replace(str, "]" => "")
                    str = replace(str, " " => "")[1:end-1]
                    strVec = split(str, ","; keepempty = false)

                    # Fill vector
                    resize!(ilbs, length(strVec))
                    ilbs .= map(s->parse(Float64, s), strVec)
                elseif occursin("Min. Neighborhood Frac:", line)
                    # Process string 
                    str = split(line, "\t"; keepempty = false)[2]
                    mnf = parse(Float64, str)
                elseif occursin("Initialization Func. Tol.:", line)
                    # Process string 
                    str = split(line, "\t"; keepempty = false)[2]
                    funcTol = parse(Float64, str) 
                elseif occursin("Initialization Max Iterations:", line)
                    str = split(line, "\t"; keepempty = false)[2]
                    maxIters = parse(Int, str)
                elseif occursin("Initialization Max Stall Iterations:", line)
                    str = split(line, "\t"; keepempty = false)[2]
                    maxSIters = Int(parse(Float64, str))
                elseif occursin("Initialization Max Time:", line)
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, " sec" => "")
                    maxTime = parse(Float64, str)
                elseif occursin("Initialization Max Stall Time:", line)
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, " sec" => "")
                    maxSTime = parse(Float64, str)
                end
            elseif inConvData
                if occursin("Time to Initialize:", line)
                    str = split(line, "\t"; keepempty = false)[2]
                    str = replace(str, "sec" => "")
                    toi = parse(Float64, str)
                elseif occursin("Initialization Iterations:", line)
                    initIters = parse(Int, split(line, "\t"; keepempty = false)[2])
                elseif occursin("Initialization Func. Evals.:", line)
                    initFEvals = parse(Int, split(line, "\t"; keepempty = false)[2])
                elseif occursin("Heuristic Opj. Function:", line)
                    hoFVal = parse(Float64, split(line, "\t"; keepempty = false)[2])
                elseif occursin("Initial Guess Converged:", line)
                    if occursin("Yes", line)
                        initConv = true
                    end
                elseif occursin("Homotopy Converged:", line)
                    if occursin("Yes", line)
                        hConv = true
                    end
                end
            elseif inCoStateData
                if occursin("# END COSTATE DATA", line)
                    inCoStateData = false
                else
                    # Read in data
                    if csDataCnt == 2
                        strVec = split(line, "\t"; keepempty = false)
                        resize!(initCS, length(strVec)) 
                        initCS .= map(s->parse(Float64, s), strVec)
                    end
                    if csDataCnt != 1 
                        if csDataCnt % 2 == 1
                            str = replace(line, "param:" => "")
                            str = replace(str, "converged:" => "")
                            strVec = split(str, " "; keepempty = false)
                            push!(hParams, parse(Float64, strVec[1]))
                            if occursin("Yes", line)
                                push!(hStepConv, true)
                            else
                                push!(hStepConv, false)
                            end
                        else
                            strVec = split(line, "\t"; keepempty = false)
                            numVec = map(s->parse(Float64, s), strVec)
                            push!(hSols, numVec)
                        end
                    end
                    
                    # Increment line counter 
                    csDataCnt += 1
                end
            end

            # Check for location in file
            if occursin("# META DATA", line)
                inMetaData = true
            elseif occursin("# CONVERGENCE DATA", line)
                inConvData = true
            elseif occursin("# COSTATE DATA", line)
                inCoStateData = true
            elseif occursin("# END META DATA", line)
                inMetaData = false
            elseif occursin("# END CONVERGENCE DATA", line) 
                inConvData = false
            end
        end

        # Push results to DataFrame
        push!(df, Dict(
            "SolutionMethod"                => solMethod,
            "InitializationCost"            => initCost,
            "CostFunctionWeights"           => cFW,
            "HeuristicOptimizer"            => ho,
            "NumberOfParticles"             => numP,
            "UsingHomotopy"                 => usingH,
            "TimeToInitialize"              => toi,
            "InitializationIterations"      => initIters,
            "InitializationFunctionEvals"   => initFEvals,
            "HeuristicObjectiveFunction"    => hoFVal,
            "InitialGuessConverged"         => initConv,
            "HomotopyConverged"             => hConv,
            "InitializedCoStates"           => initCS,
            "HomotopySolutions"             => hSols,
            "HomotopyParams"                => hParams,
            "HomotopyStepConverged"         => hStepConv,
            "InitializationUBs"             => ubs,
            "InitializationLBs"             => lbs,
            "InitializationiUBs"            => iubs,
            "InitializationiLBs"            => ilbs,
            "MinNeighborhoodFraction"       => mnf,
            "InitializationFuncTol"         => funcTol,
            "InitializationMaxIters"        => maxIters,
            "InitializationMaxStallIters"   => maxSIters,
            "InitializationMaxTime"         => maxTime,
            "InitializationMaxStallTime"    => maxSTime
        ))

        # Close files
        close(f)
    end

    return df
end

function getTextFileVersion(files)
    ver = v"0.1"
    for i in 1:length(files)
        # Open file 
        f = open(files[i], "r")

        inMetaData = false
        for line in readlines(f)

            # If in META data section, serch for version number
            if inMetaData
                if occursin("File Format Version:", line)
                    localVer = VersionNumber(split(line, "\t", keepempty = false)[2])
                    if i == 1
                        ver = localVer
                    else
                        if ver != localVer
                            throw(ErrorException("Data folder appears to contain text files with different version numbers."))
                        end
                    end
                end

                # Found version number! Break from for loop.
                break
            end

            # Search for entry and exit from META data section
            if occursin("# META DATA", line)
                inMetaData = true
            elseif ocursin("# END META DATA")
                throw(ErrorException("Looped through meta data without finding version number. " *
                    "Something is likely wrong with text file: " * files[i]))
            end
        end

        # Close file
        close(f)
    end

    return ver
end

# This method should used if encontering issues loading binary data.
# Is not able to reconstruct heuristic optimizer full solution
# but can at least grab initialized co-states.
#
# Assumes data folder contains solutions which use all similar settings
function resurrectBinaryData(folder::String, scenario::String, iConds::AbstractVector, fConds::AbstractVector,
    tspan::Tuple; homotopy = true)

    # Read in text data
    df = readTextData(folder)

    # Reconstruct indirect optimization problem
    iop = IndirectOptimizationProblem(scenario, iConds, fConds, tspan; homotopy = homotopy)

    # Reconstruct Indirect Trajectory Optimizers
    @warn "MS_PSO number of particles and number of swarms are hard coded as 500 and 8 respectively."
    dataVec = [IndirectTrajOptimizer(iop,
        solutionMethod = df.SolutionMethod[i],
        initCostFunc = df.InitializationCost[i],
        initOptimizer = df.HeuristicOptimizer[i],
        numParticles = (df.HeuristicOptimizer[i] == :PSO ? df.NumberOfParticles[i] : 500),
        numSwarms = (df.HeuristicOptimizer[i] == :PSO ? 0 : 8),
        weights = df.CostFunctionWeights[i],
        homotopyParamVec = df.HomotopyParams[i],
        dataFolder = folder * "_resurrected",
        writeData = true) for i in 1:nrow(df)]

    # Set initialized costates
    for i in 1:length(dataVec)
        initializeData!(dataVec[i].solver, df.InitializedCoStates[i])
    end
    
    return dataVec
end