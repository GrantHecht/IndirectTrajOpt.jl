using IndirectTrajOpt
using Heuristics
using IndirectCoStateInit
using IndirectShooting
using StaticArrays

# Initialize BVP function w/o STM
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
bvpFunc(y0) = cr3bpOptIntegrate(y0, tspan, ps, copyParams = true, termCallbacks = true)

# Initialize BVP function w/ STM 
psSTM = initCR3BPIndirectWithSTMParams("Low Thrust 10 CR3BP")
bvpFuncNoSTM(z0) = cr3bpOptIntegrate(z0, tspan, ps, SingleOutput(); copyParams = true)
bvpFuncWSTM(z0, ϵ) = cr3bpOptWithSTMIntegrate(z0, tspan, ϵ, psSTM; copyParams = true)

# Set initial and final conditions 
ics = [-0.019488511458668, -0.016033479812051, 0.0,
        8.918881923678198, -4.081793688818725, 0.0,
        1.0] 
fcs = [ 0.823385182067467, 0.0, -0.022277556273235,
        0.0, 0.134184170262437, 0.0, 0.0]

# Initialize FFS Initializer 
csInitializer = FSSCoStateInitializer(bvpFunc, ics, fcs; 
                                    optimizer = :MS_PSO,
                                    numParticles = 1000,
                                    initMethod = :Uniform,
                                    displayInterval = 5,
                                    maxStallIters = 50,
                                    #UBs = [40, 40, 40, 2, 2, 2, 2],
                                    #LBs = [-40, -40, -40, -2, -2, -2, -2],
                                    UBs = [100, 100, 100, 10, 10, 10, 10],
                                    LBs = [-100, -100, -100, -10, -10, -10, -10],
                                    iUBs = [40, 40, 1, 2, 2, 2, 2],
                                    iLBs = [-40, -40, -1, -2, -2, -2, -2]
                                    );

# Initialize co-states 
initialize!(csInitializer)

# Initialize Shooting Solver 
solver = FSSSolver(GetInitializedCostates(csInitializer), ics, fcs, bvpFuncNoSTM, bvpFuncWSTM; 
    homotopy = true, homotopyParamVec = [(j^2 - 1)/(25^2 - 1) for j in 25:-1:1])

solve!(solver)
