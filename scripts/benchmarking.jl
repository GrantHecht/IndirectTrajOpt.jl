using IndirectTrajOpt
using BenchmarkTools
using StaticArrays
using DifferentialEquations

# Initialize parameters
ps = initCR3BPIndirectWithSTMParams("Low Thrust CR3BP")
ps.ϵ = 0.0

# Initial condition vector
y0 = [-0.019488511458668, -0.016033479812051, 0.0,
       8.918881923678198, -4.081793688818725, 0.0, 1.0,
       15.616017, 32.875896, -0.094522,
       -0.101606, 0.044791, -0.000150, 0.133266]
y0static = SVector{14}(y0)

# Initial condition vector with STM 
z0 = [y0; zeros(196)]
for i in 1:14
    z0[14 + i + (i-1)*14] = 1.0
end

# Indirect EOMs without STM
@benchmark IndirectTrajOpt.cr3bpEomIndirect($y0, $ps, 0.0)
dy0 = zeros(14)
@benchmark IndirectTrajOpt.cr3bpEomIndirect!($dy0, $y0, $ps, 0.0)

# Indirect EOMs with STM
dz0 = zeros(210)
@benchmark IndirectTrajOpt.cr3bpEomIndirectWithSTM!($dz0, $z0, $ps, 0.0)

# Integration benchmarking
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)

# Without STM 
@benchmark cr3bpOptIntegrate($y0static, $tspan, $ps)
@benchmark cr3bpOptIntegrate($y0static, $tspan, $ps; termCallbacks = true)
@benchmark cr3bpOptIntegrate($y0, $tspan, $ps; inPlace = true)
@benchmark cr3bpOptIntegrate($y0, $tspan, $ps; inPlace = true, termCallbacks = true)

# With STM
@benchmark cr3bpOptWithSTMIntegrate($z0, $tspan, $ps)
prob = IndirectTrajOpt.createCR3BPODEWithSTMProb(z0, tspan, ps)
@benchmark solve($prob, Vern9(); reltol = 1e-14, abstol = 1e-14, save_everystep = false,
    save_start = false, initialize_save = false, maxiters = 1e6) 

