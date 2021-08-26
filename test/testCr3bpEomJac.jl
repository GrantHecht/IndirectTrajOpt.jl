
using Test
using ForwardDiff
using IndirectTrajOpt

# Initialize parameters
ps = initCR3BPIndirectWithSTMParams("Low Thrust CR3BP")

# Initial condition vector
y0 = [-0.019488511458668, -0.016033479812051, 0.0,
       8.918881923678198, -4.081793688818725, 0.0, 1.0,
       15.616017, 32.875896, -0.094522,
       -0.101606, 0.044791, -0.000150, 0.133266]

# Initial condition vector with STM 
z0 = [y0; zeros(196)]
for i in 1:14
    z0[14 + i + (i-1)*14] = 1.0
end

ps.ϵ = 1.0
for utype in [0,1,2]
    ps.utype = utype
    # Evaluate eoms with Jacobian
    dz0 = zeros(210)
    IndirectTrajOpt.cr3bpEomIndirectWithSTM!(dz0, z0, ps, 0.0)

    # Evaluate Jacobian of eoms with ForwardDiff
    dy0 = zeros(14)
    jac = ForwardDiff.jacobian(
        (y,x)->IndirectTrajOpt.cr3bpEomIndirect!(y,x,ps,0.0), 
        dy0, y0)
    ps.utype
    # Evaluate jacobian difference
    jacDiff = abs.(reshape(dz0[15:end], (14,14)) .- jac)

    # Tests 
    for i in 1:196
        @test jacDiff[i] < 1e-9
    end
end

ps.ϵ = 0.0
for utype in [0,2]
    ps.utype = utype
    # Evaluate eoms with Jacobian
    dz0 = zeros(210)
    IndirectTrajOpt.cr3bpEomIndirectWithSTM!(dz0, z0, ps, 0.0)

    # Evaluate Jacobian of eoms with ForwardDiff
    dy0 = zeros(14)
    jac = ForwardDiff.jacobian(
        (y,x)->IndirectTrajOpt.cr3bpEomIndirect!(y,x,ps,0.0), 
        dy0, y0)
    ps.utype
    # Evaluate jacobian difference
    jacDiff = abs.(reshape(dz0[15:end], (14,14)) .- jac)

    # Tests 
    for i in 1:196
        @test jacDiff[i] < 1e-9
    end
end

# Integrate state and co-states with STM 
ps.ϵ = 0.0
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
zf = cr3bpOptWithSTMIntegrate(z0, tspan, ps)

# Compute STM with ForwardDiff
stm = ForwardDiff.jacobian(
    (x)->cr3bpOptIntegrate(x, tspan, ps), y0)

# Evaluate STM difference
jacDiff = abs.(reshape(zf[15:end], (14,14)) .- stm)

# Tests 
for i in 1:196
    @test jacDiff[i] < 1e-2
end
