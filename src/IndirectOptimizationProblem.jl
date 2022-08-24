
abstract type AbstractIndirectOptimizationProblem end

# Not implementing this for now
# May want to allow users to only suppply bvps without STM and use AD for STM computation
struct UserProvidedIndirectOptimizationProblem <: AbstractIndirectOptimizationProblem
end

struct IndirectOptimizationProblem{IBVPF, BVPF, BVPSTMF} <: AbstractIndirectOptimizationProblem
    # Time span
    tspan::Tuple{Float64,Float64}

    # BVP for initialization phase (likely to include additional callbacks for early termination)
    initBVPFunc::IBVPF

    # BVP without STM 
    BVPFunc::BVPF

    # BVP with STM 
    BVPWithSTMFunc::BVPSTMF

    # Initial conditions 
    iConds::Vector{Float64}

    # Final conditions 
    fConds::Vector{Float64}

    # Flag to indicate if homotopy continuation will be employed when solving BVP. 
    # If true, BVP functions must accept continuation parameter
    homotopy::Bool
end

# TSPAN in units of days!!! Scaled for respective problems!
function IndirectOptimizationProblem(scenario::String, iConds::AbstractVector, fConds::AbstractVector, tspan::Tuple; 
    homotopy = true, homotopyFlag::HomotopyFlag = MEMF(), initializationFlag::InitializationFlag = Initialization())

    # Check for scenario type
    if occursin("CR3BP", scenario) # CR3BP Scenario
        # Initialize parameters 
        ps = initCR3BPIndirectParams(scenario)
        psSTM = initCR3BPIndirectWithSTMParams(scenario)

        # Scale tspan
        daysToTU = 86400/ps.crp.TU
        tspanScalled = (tspan[1]*daysToTU, tspan[2]*daysToTU)

        # Initialize functions
        initBVPFunc(y0, tspan) = integrate(y0, tspan, ps, CR3BP(), initializationFlag, homotopyFlag; 
                                    copyParams = true, termCallbacks = true)
        if homotopy == true
            BVPFunc(y0, tspan, ϵ)           = integrateWithHomotopy(y0, tspan, ϵ, ps, 
                                                CR3BP(), Solving(), homotopyFlag; copyParams = true)
            BVPWithSTMFunc(z0, tspan, ϵ)    = integrateWithHomotopy(z0, tspan, ϵ, psSTM, 
                                                CR3BP(), SolvingWithSTM(), homotopyFlag; copyParams = true)
        else
            BVPFunc(y0, tspan)          = integrate(y0, tspan, ps, CR3BP(), Solving(), homotopyFlag; copyParams = true)
            BVPWithSTMFunc(z0, tspan)   = integrate(z0, tspan, psSTM, CR3BP(), SolvingWithSTM(), homotopyFlag; copyParams = true)
        end
    else
        throw(ArgumentError("Only CR3BP scenarios are implemented now."))
    end

    IndirectOptimizationProblem{typeof(initBVPFunc), typeof(BVPFunc), typeof(BVPWithSTMFunc)}(
        tspanScalled, initBVPFunc, BVPFunc, BVPWithSTMFunc, iConds, fConds, homotopy
    )
end
