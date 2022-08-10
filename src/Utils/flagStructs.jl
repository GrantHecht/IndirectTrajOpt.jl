
# Dynamics model flags
abstract type DynamicsFlag end
struct CR3BP <: DynamicsFlag end

# Homotopy flags
abstract type HomotopyFlag end
struct MEMF <: HomotopyFlag end     # Min. Energy -> Min. Fuel 
struct HypTanMF <: HomotopyFlag end # Hyperbolic Tangent Minimum Fuel 

# Desired integration and returned solution
# Flag indicating integration to be performed for initialization cost function
abstract type IntegrationFlag end
# Returns full state vector and time to final desired time (for early termination)
# as first and second return values respectively.
abstract type InitializationFlag <: IntegrationFlag end
struct Initialization <: InitializationFlag end
struct InitializationWithJacobianRankPenalty <: InitializationFlag end
struct InitializationWithIntegralCost <: InitializationFlag end

# Flag indicating integration to be performed for solving BVP via 
# shooting.
#
# Returns final full state vector
abstract type SolvingFlag <: IntegrationFlag end
struct Solving <: SolvingFlag end
struct SolvingWithSTM <: SolvingFlag end

# Flag indicating full state history is disired
#
# Returns DifferentialEquations.jl Solution object
abstract type SolutionHistoryFlag <: IntegrationFlag end
struct FullSolutionHistory <: SolutionHistoryFlag end
struct FullSolutionHistoryWithSTM <: SolutionHistoryFlag end

# Flag indicating full state history without 
# control is desired
struct FullSolutionHistoryNoControl <: SolutionHistoryFlag end
# struct FullSolutionHistoryNoControlWithSTM end # Do not have this implemented!