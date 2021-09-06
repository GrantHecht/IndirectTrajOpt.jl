module IndirectTrajOpt

using IndirectCoStateInit
using IndirectShooting
import IndirectCoStateInit: initialize!
import IndirectShooting: solve!
using JLD2

# This stuff should be moved to a different package eventually
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using LoopVectorization
using Octavian
#using ModelingToolkit

# Utils
include("Utils/Spacecraft.jl")
include("Utils/matVecMulUtils.jl")
include("Utils/flagStructs.jl")

# CR3BP 
include("CR3BP/cr3bpEoms.jl")
include("CR3BP/cr3bpOptEoms.jl")
include("CR3BP/cr3bpOptEomsWithSTM.jl")
include("CR3BP/cr3bpOptEomsDiffEqCallbacks.jl")
include("CR3BP/cr3bpDiffEqUtils.jl")
include("CR3BP/cr3bpOptIntegrate.jl")
include("CR3BP/initCR3BPIndirectParams.jl")
#include("CR3BP/cr3bpEomsMTK.jl")

# Indirect Optimization 
include("IndirectOptimizationProblem.jl")
include("DataOutputManager.jl")
include("IndirectTrajOptimizer.jl")

# Temporary include. !!! Should be removed when DifferentialEquations.jl update is released !!!
#include("Utils/tempBackupIntegration.jl")

# Exports 
export SingleOutput
export initCR3BPIndirectParams
export initCR3BPIndirectWithSTMParams
export cr3bpOptIntegrate
export cr3bpOptWithSTMIntegrate
export IndirectOptimizationProblem
export IndirectTrajOptimizer
export initialize!
export solve!

end
