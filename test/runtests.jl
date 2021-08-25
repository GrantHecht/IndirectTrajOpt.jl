using IndirectTrajOpt
using Test
using SafeTestsets

@time begin
@time @safetestset "CR3BP Jacobian" begin include("testCr3bpEomJac.jl") end
end