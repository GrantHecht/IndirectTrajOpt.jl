
function cr3bpOptIntegrate(y0, tspan, ps::CR3BPIndirectParams)

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

    # Callbacks 
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
    prob = ODEProblem(ff, y0, tspan, ps; callback=cb)

    # Solve ode 
    sol = solve(
        prob,
        Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )
end