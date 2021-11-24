module IndirectTrajOpt

#using Reexport
#@reexport using IndirectCoStateInit
#@reexport using IndirectShooting
using IndirectCoStateInit
using IndirectShooting
import IndirectCoStateInit: initialize!
import IndirectShooting: solve!
using JLD2
using BSON
using DataFrames
using ProgressMeter
using Suppressor

# This stuff should be moved to a different package eventually
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using LoopVectorization
using Octavian

# Utils
include("Utils/Spacecraft.jl")
include("Utils/matVecMulUtils.jl")
include("Utils/flagStructs.jl")
include("Utils/readData.jl")

# CR3BP 
include("CR3BP/cr3bpEoms.jl")
include("CR3BP/cr3bpOptEoms.jl")
include("CR3BP/cr3bpOptEomsWithSTM.jl")
include("CR3BP/cr3bpOptEomsDiffEqCallbacks.jl")
include("CR3BP/cr3bpDiffEqUtils.jl")
include("CR3BP/cr3bpOptIntegrate.jl")
include("CR3BP/initCR3BPIndirectParams.jl")

# Indirect Optimization 
include("IndirectOptimizationProblem.jl")
include("Utils/DataOutputManager.jl")
include("IndirectTrajOptimizer.jl")
include("Utils/writeData.jl")

# Exports 
# Utility functions
export readBinaryData
export readTextData
export resurrectBinaryData

# Integration Flags
export CR3BP
export MEMF
export HypTanMF
export Initialization
export InitializationWithIntegralCost
export Solving
export SolvingWithSTM
export FullSolutionHistory
export FullSolutionHistoryWithSTM
export FullSolutionHistoryNoControl

# Parmeter initialization functions (Should be generalized with flags)
export initCR3BPIndirectParams
export initCR3BPIndirectWithSTMParams

# Integration functions
export integrate
export integrateWithHomotopy

# Indirect trajectory optimization
export IndirectOptimizationProblem
export IndirectTrajOptimizer
export initialize!
export solve!
export tSolve!

end
