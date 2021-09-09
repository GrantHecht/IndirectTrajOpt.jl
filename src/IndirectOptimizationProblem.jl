
abstract type AbstractIndirectOptimizationProblem end

# Not implementing this for now
# May want to allow users to only suppply bvps without STM and use AD for STM computation
struct UserProvidedIndirectOptimizationProblem <: AbstractIndirectOptimizationProblem
end

struct IndirectOptimizationProblem{IBVPF, BVPF, BVPSTMF} <: AbstractIndirectOptimizationProblem
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
function IndirectOptimizationProblem(scenario::String, iConds::AbstractVector, fConds::AbstractVector, 
                                     tspan::Tuple; homotopy = true)

    # Check for scenario type
    if occursin("CR3BP", scenario)
        # Initialize parameters 
        ps = initCR3BPIndirectParams(scenario)
        psSTM = initCR3BPIndirectWithSTMParams(scenario)

        # Scale tspan
        daysToTU = 86400/ps.crp.TU
        tspanScalled = (tspan[1]*daysToTU, tspan[2]*daysToTU)

        # Initialize functions
        initBVPFunc(y0) = cr3bpOptIntegrate(y0, tspanScalled, ps; copyParams = true, termCallbacks = true)
        if homotopy == true
            BVPFunc(y0, ϵ)        = cr3bpOptIntegrate(y0, tspanScalled, ϵ, ps, SingleOutput(); copyParams = true)
            BVPWithSTMFunc(z0, ϵ) = cr3bpOptWithSTMIntegrate(z0, tspanScalled, ϵ, psSTM; copyParams = true)
        else
            BVPFunc(y0)         = cr3bpOptIntegrate(y0, tspanScalled, ps, SingleOutput(); copyParams = true)
            BVPWithSTMFunc(z0)  = cr3bpOptWithSTMIntegrate(z0, tspanScalled, psSTM; copyParams = true)
        end
    else
        throw(ArgumentError("Only CR3BP scenarios are implemented now."))
    end

    IndirectOptimizationProblem{typeof(initBVPFunc), typeof(BVPFunc), typeof(BVPWithSTMFunc)}(
        initBVPFunc, BVPFunc, BVPWithSTMFunc, iConds, fConds, homotopy
    )
end