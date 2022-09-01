using IndirectTrajOpt
using Heuristics
using IndirectCoStateInit
using IndirectShooting
using StaticArrays

# Initialize BVP function w/o STM
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
bvpFunc(y0, tspan)          = integrate(y0, tspan, ps, CR3BP(), Initialization(), MEMF();
                                copyParams = true, termCallbacks = true)

# Initialize BVP function w/ STM 
psSTM = initCR3BPIndirectWithSTMParams("Low Thrust 10 CR3BP")
bvpFuncNoSTM(z0, tspan, ϵ)  = integrateWithHomotopy(z0, tspan, ϵ, ps, CR3BP(), Solving(), MEMF(); 
                                copyParams = true, inPlace = true)
bvpFuncWSTM(z0, tspan, ϵ)   = integrateWithHomotopy(z0, tspan, ϵ, psSTM, CR3BP(), SolvingWithSTM(), MEMF(); 
                                copyParams = true)

# Set initial and final conditions 
ics = [-0.019488511458668, -0.016033479812051, 0.0,
        8.918881923678198, -4.081793688818725, 0.0,
        1.0] 
fcs = [ 0.823385182067467, 0.0, -0.022277556273235,
        0.0, 0.134184170262437, 0.0, 0.0]

# Initialize FFS Initializer 
csInitializer = FSSCoStateInitializer(bvpFunc, tspan, ics, fcs; 
                                      optimizer         = :PSO,
                                      numParticles      = 500,
                                      initMethod        = :Uniform,
                                      displayInterval   = 5,
                                      maxStallIters     = 50,
                                      UBs               = [100., 100., 100., 10., 10., 10., 10.],
                                      LBs               = [-100., -100., -100., -10., -10., -10., -10.],
                                      iUBs              = [40., 40., 1., 2., 2., 2., 2.],
                                      iLBs              = [-40., -40., -1., -2., -2., -2., -2.]
                                    );

# Initialize co-states 
#initialize!(csInitializer)

# Initialize Shooting Solver 
solver = FMSSolver(tspan, ics, fcs, bvpFuncNoSTM, bvpFuncWSTM; 
    homotopy = true, homotopyParamVec = [(j^2 - 1)/(25^2 - 1) for j in 25:-1:1],
    nSeg = 2)
#solver = FSSSolver(zeros(7), tspan, ics, fcs, bvpFuncNoSTM, bvpFuncWSTM;
#    homotopy = true, homotopyParamVec = [(j^2 - 1)/(25^2 - 1) for j in 25:-1:1])

#initializeData!(solver, GetInitializedCostates(csInitializer))
initializeData!(solver, [38.001552852205165,-8.852425666012932,8.0455841251262,-0.03299197673388485,-0.1158935565431321,0.01232378525555284,3.6106167249739882])

solve!(solver)
