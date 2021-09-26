# createCR3BPODEProb : Creates a DifferentialEquations.jl ODE Problem for CR3BP state and 
# costate differential equations. Uses vector continuous callback for switch detection
# and termination if fuel is depleated for impact with body. A copy of parameters is 
# created if copyParams = true, which is required if numerically integrating in parallel
function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, params::AbstractCR3BPIndirectParams, 
    integrationTypeFlag::IntegrationFlag, homotopyFlag::MEMF; copyParams = false, 
    termCallbacks = false, inPlace = false, ϵ = -1.0, save_positions = (false, false))

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
            interp_points = 20,
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
            interp_points = 20,
            abstol = 1e-14,
            reltol = 0.0,
            save_positions = save_positions)
    end

    # ODE Problem
    if inPlace
        ff = ODEFunction{true}((du,u,p,t) -> cr3bpEomIndirect!(du,u,p,t,homotopyFlag))
    else
        ff = ODEFunction{false}((u,p,t) -> cr3bpEomIndirect(u,p,t,homotopyFlag))
    end
    return ODEProblem(ff, y0, tspan, ps; callback=cb)
end

function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, params::AbstractCR3BPIndirectParams, 
    integrationTypeFlag::InitializationWithIntegralCost, homotopyFlag::MEMF; copyParams = false, 
    termCallbacks = false, inPlace = false, ϵ = -1.0, save_positions = (false, false))

    if length(y0) != 15
        throw(ArgumentError("When computing the integral cost, y0 must be of length 15."))
    end

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
            interp_points = 20,
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
            interp_points = 20,
            abstol = 1e-14,
            reltol = 0.0,
            save_positions = save_positions)
    end

    # ODE Problem
    if inPlace
        ff = ODEFunction{true}((du,u,p,t) -> cr3bpEomIndirectIntegralCost!(du,u,p,t,homotopyFlag))
    else
        ff = ODEFunction{false}((u,p,t) -> cr3bpEomIndirectIntegralCost(u,p,t,homotopyFlag))
    end
    return ODEProblem(ff, y0, tspan, ps; callback=cb)
end


function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, scenario::String, 
    integrationTypeFlag::IntegrationFlag, homotopyFlag::HomotopyFlag)
    ps = initCR3BPIndirectParams(scenario)
    return createCR3BPODEProb(y0, tspan, ps, integrationTypeFlag, homotopyFlag)
end

function createCR3BPODEWithSTMProb(z0::AbstractVector, tspan::Tuple, params::CR3BPIndirectWithSTMParams, 
    homotopyFlag::MEMF; copyParams = false, ϵ = -1.0)

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
    ff = ODEFunction{true}((du,u,p,t) -> cr3bpEomIndirectWithSTM!(du,u,p,t,homotopyFlag))
    return ODEProblem(ff, z0, tspan, ps; callback=cb)
end

function createCR3BPODEWithSTMProb(z0::AbstractVector, tspan::Tuple, scenario::String, homotopyFlag::HomotopyFlag)
    ps = initCR3BPIndirectWithSTMParams(scenario)
    return createCR3BPODEWithSTMProb(z0, tspan, ps, homotopyFlag)
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