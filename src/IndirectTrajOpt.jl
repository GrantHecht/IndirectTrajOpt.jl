module IndirectTrajOpt

using StaticArrays
using LinearAlgebra
using DifferentialEquations
using ModelingToolkit

# Utils
include("Spacecraft.jl")

# CR3BP 
include("CR3BP/cr3bpEoms.jl")
include("CR3BP/cr3bpOptEoms.jl")
include("CR3BP/cr3bpOptEomsDiffEqCallbacks.jl")
include("CR3BP/cr3bpOptIntegrate.jl")
include("CR3BP/initCR3BPIndirectParams.jl")
include("CR3BP/cr3bpEomsMTK.jl")

# Exports 
export initCR3BPIndirectParams
export cr3bpOptIntegrate

end
