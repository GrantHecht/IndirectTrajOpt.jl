# CR3BP Indirect with STM parameters 
mutable struct CR3BPIndirectWithSTMParams <: AbstractCR3BPIndirectParams
    # Spacecraft parameters 
    sp::SimpleSpacecraft

    # CR3BP parameters
    crp::CR3BPParams

    # Continuation variable
    ϵ::Float64

    # Thrust type 
    utype::Int64

    # Allocate matricies
    m1::Matrix{Float64}
    m2::Matrix{Float64}

    function CR3BPIndirectWithSTMParams(sp::Spacecraft, crp::CR3BPParams, ϵ)
        new(sp, crp, ϵ, 0, zeros(14,14), zeros(14,14))
    end
end

function cr3bpEomIndirectWithSTM!(du::AbstractVector, u::AbstractVector,
                                  p::CR3BPIndirectWithSTMParams, t,
                                  homotopyFlag::MEMF)
    @inbounds begin
        # Evaluate state/co-state dynamics
        GVec = cr3bpEomIndirect!(du, u, p, t, homotopyFlag)

        # Dynamics Jacobian
        jac = p.m1; jac .= 0.0
        cr3bpEOMJac!(jac, u, p, GVec, homotopyFlag)

        # Compute STM derivs 
        Φ   = reshape(view(u,  15:210), (14,14))
        dΦ  = reshape(view(du, 15:210), (14,14))
        Octavian.matmul!(dΦ, jac, Φ, 1.0, 0.0)
    end
    return nothing
end

function cr3bpEOMJac!(jac, u, p, GVec, homotopyFlag::MEMF)

    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    MU  = p.crp.MU 
    μ   = p.crp.μ
    c   = p.sp.c

    # Scale requirements
    tMaxSc = p.sp.tMax * TU * TU / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p.utype, p.ϵ)

    # Compute r1 and r2 with inverse powers
    xpmu    = u[1] + μ
    xpmum1  = u[1] + μ - 1
    r1      = sqrt(xpmu*xpmu + u[2]*u[2] + u[3]*u[3])
    r2      = sqrt(xpmum1*xpmum1 + u[2]*u[2] + u[3]*u[3])
    invr1   = inv(r1)
    invr2   = inv(r2)
    invr15  = invr1*invr1*invr1*invr1*invr1
    invr25  = invr2*invr2*invr2*invr2*invr2
    invr17  = invr15*invr1*invr1
    invr27  = invr25*invr2*invr2

    # Partials of G matrix wrt r
    dG11dr          = @SVector [(1 - μ)*(μ + u[1])*(9*invr15 - 15*(u[1] + μ)^2*invr17) +
                    μ*(μ + u[1] - 1)*(9*invr25 - 15*(u[1] + μ - 1)^2*invr27),
                    (1 - μ)*u[2]*(3*invr15 - 15*(u[1] + μ)^2*invr17) +
                    μ*u[2]*(3*invr25 - 15*(u[1] + μ - 1)^2*invr27), 
                    (1 - μ)*u[3]*(3*invr15 - 15*(u[1] + μ)^2*invr17) + 
                    μ*u[3]*(3*invr25 - 15*(u[1] + μ - 1)^2*invr27)];
                    
    dG12dr          = @SVector [u[2]*(1 - μ)*(3*invr15 - 15*(u[1] + μ)^2*invr17) + 
                    u[2]*μ*(3*invr25 - 15*(u[1] + μ - 1)^2*invr27), 
                    (1 - μ)*(u[1] + μ)*(3*invr15 - 15*u[2]^2*invr17) + 
                    μ*(u[1] + μ - 1)*(3*invr25 - 15*u[2]^2*invr27), 
                    -15*(1 - μ)*(u[1] + μ)*u[2]*u[3]*invr17 - 
                        15*μ*(u[1] + μ - 1)*u[2]*u[3]*invr27];

    dG13dr          = @SVector [u[3]*(1 - μ)*(3*invr15 - 15*(u[1] + μ)^2*invr17) + 
                    u[3]*μ*(3*invr25 - 15*(u[1] + μ - 1)^2*invr27), 
                    -15*(1 - μ)*(u[1] + μ)*u[2]*u[3]*invr17 - 
                    15*μ*(u[1] + μ - 1)*u[2]*u[3]*invr27, 
                    (1 - μ)*(u[1] + μ)*(3*invr15 - 15*u[3]^2*invr17) + 
                    μ*(u[1] + μ - 1)*(3*invr25 - 15*u[3]^2*invr27)];

    dG22dr          = @SVector [(1 - μ)*(u[1] + μ)*(3*invr15 - 15*u[2]^2*invr17) + 
                    μ*(u[1] + μ - 1)*(3*invr25 - 15*u[2]^2*invr27), 
                    (1 - μ)*u[2]*(9*invr15 - 15*u[2]^2*invr17) + 
                    μ*u[2]*(9*invr25 - 15*u[2]^2*invr27), 
                    (1 - μ)*u[3]*(3*invr15 - 15*u[2]^2*invr17) + 
                    μ*u[3]*(3*invr25 - 15*u[2]^2*invr27)];

    dG23dr          = @SVector [-15*(1 - μ)*(u[1] + μ)*u[2]*u[3]*invr17 - 
                    15*μ*(u[1] + μ - 1)*u[2]*u[3]*invr27, 
                    (1 - μ)*u[3]*(3*invr15 - 15*u[2]^2*invr17) + 
                    μ*u[3]*(3*invr25 - 15*u[2]^2*invr27), 
                    (1 - μ)*u[2]*(3*invr15 - 15*u[3]^2*invr17) + 
                    μ*u[2]*(3*invr25 - 15*u[3]^2*invr27)];

    dG33dr          = @SVector [(1 - μ)*(u[1] + μ)*(3*invr15 - 15*u[3]^2*invr17) + 
                    μ*(u[1] + μ - 1)*(3*invr25 - 15*u[3]^2*invr27), 
                    (1 - μ)*u[2]*(3*invr15 - 15*u[3]^2*invr17) + 
                    μ*u[2]*(3*invr25 - 15*u[3]^2*invr27), 
                    (1 - μ)*u[3]*(9*invr15 - 15*u[3]^2*invr17) + 
                    μ*u[3]*(9*invr25 - 15*u[3]^2*invr27)];  

    # Fill Jacobian
    jac[8,   1] = -dG11dr[1]*u[11] - dG12dr[1]*u[12] - dG13dr[1]*u[13]
    jac[8,   2] = -dG11dr[2]*u[11] - dG12dr[2]*u[12] - dG13dr[2]*u[13]
    jac[8,   3] = -dG11dr[3]*u[11] - dG12dr[3]*u[12] - dG13dr[3]*u[13]
    jac[9,   1] = -dG12dr[1]*u[11] - dG22dr[1]*u[12] - dG23dr[1]*u[13]
    jac[9,   2] = -dG12dr[2]*u[11] - dG22dr[2]*u[12] - dG23dr[2]*u[13]
    jac[9,   3] = -dG12dr[3]*u[11] - dG22dr[3]*u[12] - dG23dr[3]*u[13]
    jac[10,  1] = -dG13dr[1]*u[11] - dG23dr[1]*u[12] - dG33dr[1]*u[13]
    jac[10,  2] = -dG13dr[2]*u[11] - dG23dr[2]*u[12] - dG33dr[2]*u[13]
    jac[10,  3] = -dG13dr[3]*u[11] - dG23dr[3]*u[12] - dG33dr[3]*u[13]
    jac[1,   4] = 1.0
    jac[2,   5] = 1.0
    jac[3,   6] = 1.0
    jac[4,   1] = GVec[1]
    jac[4,   2] = GVec[4]
    jac[4,   3] = GVec[5]
    jac[5,   1] = GVec[4]
    jac[5,   2] = GVec[2]
    jac[5,   3] = GVec[6]
    jac[6,   1] = GVec[5]
    jac[6,   2] = GVec[6]
    jac[6,   3] = GVec[3] 
    jac[4,   5] = 2.0
    jac[5,   4] = -2.0
    jac[8,  11] = -GVec[1]
    jac[8,  12] = -GVec[4]
    jac[8,  13] = -GVec[5]
    jac[9,  11] = -GVec[4]
    jac[9,  12] = -GVec[2]
    jac[9,  13] = -GVec[6]
    jac[10, 11] = -GVec[5]
    jac[10, 12] = -GVec[6]
    jac[10, 13] = -GVec[3] 
    jac[11,  8] = -1.0
    jac[12,  9] = -1.0
    jac[13, 10] = -1.0
    jac[11, 12] = 2.0
    jac[12, 11] = -2.0

    # Fill remaining depending on utype
    if p.utype != 1
        compute_G1!(jac, u, λv, tMaxSc, γ)
        compute_G2!(jac, u, λv, tMaxSc, γ)
        compute_G7!(jac, u, λv, tMaxSc, γ)
        compute_G8!(jac, u, λv, tMaxSc, γ)
    else
        compute_G1!(jac, u, λv, tMaxSc, cSc, p.ϵ, γ)
        compute_G2!(jac, u, λv, tMaxSc, cSc, p.ϵ, γ)
        compute_G3!(jac, u, λv, tMaxSc, p.ϵ)
        compute_G4!(jac, u, λv, tMaxSc, p.ϵ)
        compute_G5!(jac, u, λv, tMaxSc, p.ϵ)
        compute_G6!(jac, p.ϵ, tMaxSc, cSc)
        compute_G7!(jac, u, λv, tMaxSc, cSc, p.ϵ, γ)
        compute_G8!(jac, u, λv, tMaxSc, cSc, p.ϵ, γ)
        compute_G9!(jac, u, λv, tMaxSc, p.ϵ)
    end
    return nothing
end

function compute_G1!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, cSc, ϵ, γ)
    @inbounds begin
        temp = γ*tMaxSc / (λv*u[7]*u[7]) + cSc*tMaxSc / (2.0*ϵ*u[7]^3)
        jac[4, 7] = temp*u[11]
        jac[5, 7] = temp*u[12]
        jac[6, 7] = temp*u[13]
    end
    return nothing
end

function compute_G1!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, γ)
    @inbounds begin
        temp    = γ*tMaxSc / (λv*u[7]*u[7])
        jac[4, 7] = temp*u[11]
        jac[5, 7] = temp*u[12]
        jac[6, 7] = temp*u[13]
    end
    return nothing
end

function compute_G2!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, γ)
    @inbounds begin
        λvInv    = 1.0 / λv
        λvInv3   = λvInv^3
        temp1    = -γ*tMaxSc/u[7]
        jac[4, 11] = temp1*(λvInv - u[11]*u[11]*λvInv3)
        jac[5, 12] = temp1*(λvInv - u[12]*u[12]*λvInv3)
        jac[6, 13] = temp1*(λvInv - u[13]*u[13]*λvInv3)
        jac[4, 12] = -temp1*u[11]*u[12]*λvInv3 
        jac[4, 13] = -temp1*u[11]*u[13]*λvInv3
        jac[5, 11] = jac[4, 12]
        jac[5, 13] = -temp1*u[12]*u[13]*λvInv3
        jac[6, 11] = jac[4, 13]
        jac[6, 12] = jac[5, 13]
    end
    return nothing
end

function compute_G2!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, cSc, ϵ, γ)
    @inbounds begin
        λvInv   = 1.0 / λv
        λvInv2  = λvInv*λvInv
        λvInv3  = λvInv2*λvInv
        temp1   = -γ*tMaxSc / u[7]
        temp2   = -cSc*tMaxSc / (2.0*ϵ*u[7]^2)
        jac[4, 11] = temp1*(λvInv - u[11]*u[11]*λvInv3) + temp2*u[11]*u[11]*λvInv2
        jac[5, 12] = temp1*(λvInv - u[12]*u[12]*λvInv3) + temp2*u[12]*u[12]*λvInv2
        jac[6, 13] = temp1*(λvInv - u[13]*u[13]*λvInv3) + temp2*u[13]*u[13]*λvInv2
        jac[4, 12] = -temp1*u[11]*u[12]*λvInv3 + temp2*u[11]*u[12]*λvInv2 
        jac[4, 13] = -temp1*u[11]*u[13]*λvInv3 + temp2*u[11]*u[13]*λvInv2
        jac[5, 11] = jac[4, 12]
        jac[5, 13] = -temp1*u[12]*u[13]*λvInv3 + temp2*u[12]*u[13]*λvInv2
        jac[6, 11] = jac[4, 13]
        jac[6, 12] = jac[5, 13]
    end
    return nothing
end

function compute_G3!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, ϵ)
    @inbounds begin
        temp        = -tMaxSc / (2.0*ϵ*λv*u[7])
        jac[4, 14]    = temp*u[11]
        jac[5, 14]    = temp*u[12]
        jac[6, 14]    = temp*u[13]
    end
    return nothing
end

function compute_G4!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, ϵ)
    @inbounds jac[7,7] = λv*tMaxSc / (2.0*ϵ*u[7]^2)
    return nothing
end

function compute_G5!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, ϵ)
    @inbounds begin
        temp = -tMaxSc / (2.0*ϵ*λv*u[7])
        jac[7, 11] = temp*u[11]
        jac[7, 12] = temp*u[12]
        jac[7, 13] = temp*u[13]
    end
    return nothing
end

function compute_G6!(jac::AbstractMatrix, ϵ, tMaxSc, cSc)
    @inbounds jac[7,14] = -tMaxSc / (2.0*ϵ*cSc)
    return nothing
end

function compute_G7!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, γ)
    @inbounds jac[14, 7] = 2.0*λv*γ*tMaxSc / (u[7]^3)
    return nothing
end

function compute_G7!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, cSc, ϵ, γ)
    @inbounds jac[14,7] = 2.0*λv*γ*tMaxSc / (u[7]^3) + λv^2*cSc*tMaxSc / (2.0*ϵ*u[7]^4)
    return nothing
end

function compute_G8!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, γ)
    @inbounds begin
        temp = -γ*tMaxSc / (λv*u[7]^2)      
        jac[14, 11] = temp*u[11]
        jac[14, 12] = temp*u[12]
        jac[14, 13] = temp*u[13]
    end
    return nothing
end

function compute_G8!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, cSc, ϵ, γ)
    @inbounds begin
        temp = -γ*tMaxSc / (λv*u[7]^2) - cSc*tMaxSc / (2.0*ϵ*u[7]^3)
        jac[14, 11] = temp*u[11]
        jac[14, 12] = temp*u[12]
        jac[14, 13] = temp*u[13]
    end
    return nothing
end

function compute_G9!(jac::AbstractMatrix, u::AbstractVector, λv, tMaxSc, ϵ)
    @inbounds jac[14,14] = -λv*tMaxSc / (2.0*ϵ*u[7]^2)
    return nothing
end
