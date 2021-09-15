
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using BenchmarkTools
using DifferentialEquations
using ModelingToolkit
using Plots


# Initialize parameters
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
ps.ϵ = 0.0

# Initial condition vector
y0 = @SVector [-0.019488511458668, -0.016033479812051, 0.0,
                8.918881923678198, -4.081793688818725, 0.0, 1.0,
                -24.3123620400503, -60.46766521717989, -0.2883477429418446, 0.18005551834799977, -0.09283727657419928, -0.00021162595932502253, 0.5613066482340945]

# Integrate
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
#@btime cr3bpOptIntegrate($y0, $tspan, $ps)
sol = cr3bpOptIntegrate(y0, tspan, ps)

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
cb = ContinuousCallback(
    IndirectTrajOpt.cr3bpEomsConditionNoTerm,
    IndirectTrajOpt.cr3bpEomsAffectNoTerm!,
    IndirectTrajOpt.cr3bpEomsAffectNoTerm!;
    idxs = nothing,
    rootfind = true,
    interp_points = 10,
    abstol = 1e-14,
    reltol = 0.0,
    save_positions = (true, true))

# ODE Problem
ff = ODEFunction{true}(IndirectTrajOpt.cr3bpEomIndirect!)
prob = ODEProblem(ff, Vector{Float64}(y0), tspan, ps; callback=cb)

# Solve ode 
sol = solve(
    prob,
    Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
    save_everystep = true,
    save_start = true,
    initialize_save = true,
    maxiters = 1e6
    )

display(plot(sol, vars=(1,2)))

#@btime IndirectTrajOpt.cr3bpEomIndirect($y0, $ps, 0.0)
#@btime IndirectTrajOpt.cr3bpCostateEom($y0, $ps, 1.0)
#at = @SVector [1.0, 0.0, 0.0]
#@btime IndirectTrajOpt.cr3bpEomControl($y0, $ps.crp,0.0,$at)
#dy0 = similar(y0)
#@btime IndirectTrajOpt.cr3bpEomIndirect!($dy0,$y0,$ps,0.0)
