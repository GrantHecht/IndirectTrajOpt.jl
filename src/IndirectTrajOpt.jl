module IndirectTrajOpt

using StaticArrays
using LinearAlgebra
using DifferentialEquations
using LoopVectorization
using Octavian
#using ModelingToolkit

# Utils
include("Spacecraft.jl")
include("matVecMulUtils.jl")

# CR3BP 
include("CR3BP/cr3bpEoms.jl")
include("CR3BP/cr3bpOptEoms.jl")
include("CR3BP/cr3bpOptEomsWithSTM.jl")
include("CR3BP/cr3bpOptEomsDiffEqCallbacks.jl")
include("CR3BP/cr3bpDiffEqUtils.jl")
include("CR3BP/cr3bpOptIntegrate.jl")
include("CR3BP/initCR3BPIndirectParams.jl")
#include("CR3BP/cr3bpEomsMTK.jl")

# Exports 
export initCR3BPIndirectParams
export initCR3BPIndirectWithSTMParams
export cr3bpOptIntegrate
export cr3bpOptWithSTMIntegrate

end
