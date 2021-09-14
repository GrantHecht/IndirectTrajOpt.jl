
# Dynamics model flags
struct CR3BP end

# Desired integration and returned solution
# Flag indicating integration to be performed for initialization cost function
#
# Returns full state vector and time to final desired time (for early termination)
# as first and second return values respectively.
struct Initialization end

# Flag indicating integration to be performed for solving BVP via 
# shooting.
#
# Returns final full state vector
struct Solving end
struct SolvingWithSTM end

# Flag indicating full state history is disired
#
# Returns DifferentialEquations.jl Solution object
struct FullSolutionHistory end
struct FullSolutionHistoryWithSTM end

# Flag indicating full state history without 
# control is desired
struct FullSolutionHistoryNoControl end
# struct FullSolutionHistoryNoControlWithSTM end # Do not have this implemented!