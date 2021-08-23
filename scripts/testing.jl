
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using ModelingToolkit
#using Plots

function main()

    # Initialize parameters
    ps = initCR3BPIndirectParams("Low Thrust CR3BP")
    ps.ϵ = 0.0

    # Initial condition vector
    y0 = @SVector [-0.019488511458668, -0.016033479812051, 0.0,
                    8.918881923678198, -4.081793688818725, 0.0, 1.0,
                    15.616017, 32.875896, -0.094522,
                    -0.101606, 0.044791, -0.000150, 0.133266]

    # Integrate
    tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
    #@btime cr3bpOptIntegrate($y0, $tspan, $ps)
    #sol = cr3bpOptIntegrate(y0, tspan, ps)
    #display(plot(sol, vars=(1,2)))

    # Check that utype is set appropriately
    cSc = ps.sp.isp*9.81*ps.crp.TU / (ps.crp.LU*1000.0)
    λv = norm(view(y0,11:13))
    S = IndirectTrajOpt.computeS(y0, λv, cSc)
    if S > ps.ϵ
        ps.utype = 0
    elseif S < -ps.ϵ
        ps.utype = 2
    else
        ps.utype = 1
    end

    # Callbacks 
    cb = VectorContinuousCallback(
        IndirectTrajOpt.cr3bpEomsCondition,
        IndirectTrajOpt.cr3bpEomsAffect!,
        IndirectTrajOpt.cr3bpEomsAffect!, 4;
        idxs = nothing,
        rootfind = true,
        interp_points = 10,
        abstol = 1e-14,
        reltol = 0.0,
        save_positions = (false, false))

    # ODE Problem
    ff = ODEFunction{true}(IndirectTrajOpt.cr3bpEomIndirect!)
    prob = ODEProblem(ff, Vector{Float64}(y0), tspan, ps; callback=cb)

    # Solve ode 
    @btime solve(
        $prob,
        $Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
        save_everystep = false,
        save_start = false,
        initialize_save = false,
        maxiters = 1e6
        )

    #@btime IndirectTrajOpt.cr3bpEomIndirect($y0, $ps, 0.0)
    #@btime IndirectTrajOpt.cr3bpCostateEom($y0, $ps, 1.0)
    #at = @SVector [1.0, 0.0, 0.0]
    #@btime IndirectTrajOpt.cr3bpEomControl($y0, $ps.crp,0.0,$at)
    #dy0 = similar(y0)
    #@btime IndirectTrajOpt.cr3bpEomIndirect!($dy0,$y0,$ps,0.0)

end

res = main()