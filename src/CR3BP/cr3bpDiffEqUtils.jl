# createCR3BPODEProb : Creates a DifferentialEquations.jl ODE Problem for CR3BP state and 
# costate differential equations. Uses vector continuous callback for switch detection
# and termination if fuel is depleated for impact with body. A copy of parameters is 
# created if copyParams = true, which is required if numerically integrating in parallel
function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, params::CR3BPIndirectParams; copyParams = false)

    # Copy parameters if desired
    if copyParams 
        ps = deepcopy(params)
    else
        ps = params 
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
    cb = VectorContinuousCallback(
        cr3bpEomsCondition,
        cr3bpEomsAffect!,
        cr3bpEomsAffect!, 4;
        idxs = nothing,
        rootfind = true,
        interp_points = 10,
        abstol = 1e-14,
        reltol = 0.0,
        save_positions = (false, false))

    # ODE Problem
    ff = ODEFunction{false}(cr3bpEomIndirect)
    return ODEProblem(ff, y0, tspan, ps; callback=cb)
end

function createCR3BPODEProb(y0::AbstractVector, tspan::Tuple, scenario::String)
    ps = initCR3BPIndirectParams(scenario)
    return createCR3BPODEProb(y0, tspan, ps)
end

