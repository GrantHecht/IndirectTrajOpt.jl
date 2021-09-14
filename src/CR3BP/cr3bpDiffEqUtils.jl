# createCR3BPODEProb : Creates a DifferentialEquations.jl ODE Problem for CR3BP state and 
# costate differential equations. Uses vector continuous callback for switch detection
# and termination if fuel is depleated for impact with body. A copy of parameters is 
# created if copyParams = true, which is required if numerically integrating in parallel
function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, params::AbstractCR3BPIndirectParams; 
    copyParams = false, termCallbacks = false, inPlace = false, ϵ = -1.0, save_positions = (false, false))

    # Copy parameters if desired
    if copyParams 
        ps = deepcopy(params)
    else
        ps = params 
    end

    # Set ϵ if desired 
    if ϵ >= 0
        ps.ϵ = ϵ
    end
    
    # Check that utype is set appropriately
    cSc = ps.sp.isp*9.81*ps.crp.TU / (ps.crp.LU*1000.0)
    λv = norm(view(y0,11:13))
    S = computeS(y0, λv, cSc)
    if S > ps.ϵ
        ps.utype = 0
    elseif S < -ps.ϵ
        ps.utype = 2
    else
        ps.utype = 1
    end

    # Instantiate callback
    if termCallbacks
        cb = VectorContinuousCallback(
            cr3bpEomsCondition,
            cr3bpEomsAffect!,
            cr3bpEomsAffect!, 4;
            idxs = nothing,
            rootfind = true,
            interp_points = 10,
            abstol = 1e-14,
            reltol = 0.0,
            save_positions = save_positions)
    else
        cb = ContinuousCallback(
            cr3bpEomsConditionNoTerm,
            cr3bpEomsAffectNoTerm!,
            cr3bpEomsAffectNoTerm!;
            idxs = nothing,
            rootfind = true,
            interp_points = 10,
            abstol = 1e-14,
            reltol = 0.0,
            save_positions = save_positions)
    end

    # ODE Problem
    if inPlace
        ff = ODEFunction{true}(cr3bpEomIndirect!)
    else
        ff = ODEFunction{false}(cr3bpEomIndirect)
    end
    return ODEProblem(ff, y0, tspan, ps; callback=cb)
end

function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, scenario::String)
    ps = initCR3BPIndirectParams(scenario)
    return createCR3BPODEProb(y0, tspan, ps)
end

function createCR3BPODEWithSTMProb(z0::AbstractVector, tspan::Tuple, params::CR3BPIndirectWithSTMParams; 
    copyParams = false, ϵ = -1.0)

    # Copy parameters if desired
    if copyParams 
        ps = deepcopy(params)
    else
        ps = params 
    end

    # Set ϵ if desired 
    if ϵ >= 0.0
        ps.ϵ = ϵ
    end
    
    # Check that utype is set appropriately
    cSc = ps.sp.isp*9.81*ps.crp.TU / (ps.crp.LU*1000.0)
    λv = norm(view(z0,11:13))
    S = computeS(z0, λv, cSc)
    if S > ps.ϵ
        ps.utype = 0
    elseif S < -ps.ϵ
        ps.utype = 2
    else
        ps.utype = 1
    end

    # Instantiate callback
    cb = ContinuousCallback(
        cr3bpEomsConditionNoTerm,
        cr3bpEomsAffectNoTermWithSTM!,
        cr3bpEomsAffectNoTermWithSTM!;
        idxs = nothing,
        rootfind = true,
        interp_points = 10,
        abstol = 1e-14,
        reltol = 0.0,
        save_positions = (false, false))

    # ODE Problem
    ff = ODEFunction{true}(cr3bpEomIndirectWithSTM!)
    return ODEProblem(ff, z0, tspan, ps; callback=cb)
end

function createCR3BPODEWithSTMProb(z0::AbstractVector, tspan::Tuple, scenario::String)
    ps = initCR3BPIndirectWithSTMParams(scenario)
    return createCR3BPODEWithSTMProb(z0, tspan, ps)
end

function createCR3BPODENoControlProb(x0::AbstractVector, tspan::Tuple, params::CR3BPParams; 
    copyParams = false, inPlace = false)

    # Copy parameters if desired
    if copyParams 
        ps = deepcopy(params)
    else
        ps = params 
    end

   # ODE Problem
    if inPlace
        ff = ODEFunction{true}(cr3bpEomNoControl!)
    else
        ff = ODEFunction{false}(cr3bpEomNoControl)
    end

    return ODEProblem(ff, x0, tspan, ps)
end

function createCR3BPODENoControlProb(x0::AbstractVector, tspan::Tuple, scenario::String; inPlace = false)
    psFull = initCR3BPIndirectParams(scenario)
    return createCR3BPODENoControlProb(x0, tspan, psFull.crp; copyParams = false, inPlace = inPlace)
end