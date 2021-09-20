abstract type AbstractCR3BPIndirectParams end

# CR3BP Indirect parameters 
mutable struct CR3BPIndirectParams <: AbstractCR3BPIndirectParams
    # Spacecraft parameters 
    sp::SimpleSpacecraft

    # CR3BP parameters
    crp::CR3BPParams

    # Continuation variable
    ϵ::Float64

    # Thrust type 
    utype::Int64

    function CR3BPIndirectParams(sp::Spacecraft, crp::CR3BPParams, ϵ)
        new(sp, crp, ϵ, 0)
    end
end

# CR3BP CoState Dynamics
function cr3bpCostateEom(u::AbstractArray, p::AbstractCR3BPIndirectParams, γ)
    @inbounds begin
        # Get requirements 
        TU  = p.crp.TU
        LU  = p.crp.LU
        MU  = p.crp.MU 
        μ   = p.crp.μ

        # Scale requirements
        tMaxSc = p.sp.tMax * TU * TU / (MU*LU*1000.0)

        # Compute Requirements
        xpmu    = u[1] + μ
        xpmum1  = u[1] + μ - 1
        r1      = sqrt(xpmu*xpmu + u[2]*u[2] + u[3]*u[3])
        r2      = sqrt(xpmum1*xpmum1 + u[2]*u[2] + u[3]*u[3])
        invr13  = r1^(-3)
        invr15  = r1^(-5)
        invr23  = r2^(-3)
        invr25  = r2^(-5)
        λv      = sqrt(u[11]^2 + u[12]^2 + u[13]^2)

        # Compute G 
        G11     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*(u[1] + μ)^2*invr15 - 
                    μ*invr23 + 3*μ*(u[1] + μ - 1)^2*invr25;
        G22     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*u[2]^2*invr15 - 
                    μ*invr23 + 3*μ*u[2]^2*invr25;
        G33     = -(1 - μ)*invr13 + 3*(1 - μ)*u[3]^2*invr15 - 
                    μ*invr23 + 3*μ*u[3]^2*invr25;
        G12     = 3*(1 - μ)*(u[1] + μ)*u[2]*invr15 + 
                    3*μ*(u[1] + μ - 1)*u[2]*invr25;
        G13     = 3*(1 - μ)*(u[1] + μ)*u[3]*invr15 + 
                    3*μ*(u[1] + μ - 1)*u[3]*invr25;
        G23     = 3*(1 - μ)*u[2]*u[3]*invr15 + 
                    3*μ*u[2]*u[3]*invr25;

        # Compute and return Dynamics
        dλ = @SVector [ -G11*u[11] - G12*u[12] - G13*u[13],
                        -G12*u[11] - G22*u[12] - G23*u[13],
                        -G13*u[11] - G23*u[12] - G33*u[13],
                        -u[8] + 2.0*u[12],
                        -u[9] - 2.0*u[11],
                        -u[10],
                        -λv*γ*tMaxSc / (u[7]*u[7])]
    end

    return dλ
end

function cr3bpCostateEom!(dλ::AbstractArray, u::AbstractArray, 
                          p::AbstractCR3BPIndirectParams, γ)
    @inbounds begin
    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    MU  = p.crp.MU 
    μ   = p.crp.μ

    # Scale requirements
    tMaxSc = p.sp.tMax * TU^2 / (MU*LU*1000.0)

    # Compute Requirements
    r1      = sqrt((u[1] + μ)^2 + u[2]^2 + u[3]^2)
    r2      = sqrt((u[1] + μ - 1)^2 + u[2]^2 + u[3]^2)
    invr13  = r1^(-3)
    invr15  = r1^(-5)
    invr23  = r2^(-3)
    invr25  = r2^(-5)
    λv      = sqrt(u[11]^2 + u[12]^2 + u[13]^2)

    # Compute G 
    G11     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*(u[1] + μ)^2*invr15 - 
                μ*invr23 + 3*μ*(u[1] + μ - 1)^2*invr25;
    G22     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*u[2]^2*invr15 - 
                μ*invr23 + 3*μ*u[2]^2*invr25;
    G33     = -(1 - μ)*invr13 + 3*(1 - μ)*u[3]^2*invr15 - 
                μ*invr23 + 3*μ*u[3]^2*invr25;
    G12     = 3*(1 - μ)*(u[1] + μ)*u[2]*invr15 + 
                3*μ*(u[1] + μ - 1)*u[2]*invr25;
    G13     = 3*(1 - μ)*(u[1] + μ)*u[3]*invr15 + 
                3*μ*(u[1] + μ - 1)*u[3]*invr25;
    G23     = 3*(1 - μ)*u[2]*u[3]*invr15 + 
                3*μ*u[2]*u[3]*invr25;

    # Compute and return Dynamics
    dλ[1] = -G11*u[11] - G12*u[12] - G13*u[13]
    dλ[2] = -G12*u[11] - G22*u[12] - G23*u[13]
    dλ[3] = -G13*u[11] - G23*u[12] - G33*u[13]
    dλ[4] = -u[8] + 2.0*u[12]
    dλ[5] = -u[9] - 2.0*u[11]
    dλ[6] = -u[10]
    dλ[7] = -λv*γ*tMaxSc / (u[7]^2)
    end

    return @SVector [G11, G22, G33, G12, G13, G23]
end

# CR3BP Indirect EOMs 
function cr3bpEomIndirect(u::AbstractVector,p::AbstractCR3BPIndirectParams,t)
    @inbounds begin

    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    VU  = p.crp.VU 
    MU  = p.crp.MU 
    isp = p.sp.isp 
    c   = p.sp.c

    # Scale requirements
    tMaxSc = p.sp.tMax * TU * TU / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p.utype, p.ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Derivatives
    dx = cr3bpEomControl(u,p.crp,t,at)
    dλ = cr3bpCostateEom(u,p, γ)
    du = @SVector [dx[1], dx[2], dx[3], dx[4], dx[5], dx[6], -γ*tMaxSc / cSc,
                   dλ[1], dλ[2], dλ[3], dλ[4], dλ[5], dλ[6], dλ[7]]

    return du
    end
end

function cr3bpEomIndirect!(du::AbstractVector, u::AbstractVector,
                           p::AbstractCR3BPIndirectParams, t)
    @inbounds begin
    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    VU  = p.crp.VU 
    MU  = p.crp.MU 
    isp = p.sp.isp 
    c   = p.sp.c

    # Scale requirements
    tMaxSc = p.sp.tMax * TU^2 / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]^2 + u[12]^2 + u[13]^2)
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p.utype, p.ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Compute dynamics
    cr3bpEomControl!(view(du,1:6), u, p.crp, t, at)
    du[7] = -γ*tMaxSc / cSc
    GVec  = cr3bpCostateEom!(view(du,8:14), u, p, γ)
    end
    return GVec
end

function cr3bpEomIndirectEnergyIntegral(u::AbstractVector, p::AbstractCR3BPIndirectParams, t)
    @inbounds begin

    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    VU  = p.crp.VU 
    MU  = p.crp.MU 
    isp = p.sp.isp 
    c   = p.sp.c

    # Scale requirements
    tMaxSc = p.sp.tMax * TU * TU / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p.utype, p.ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Derivatives
    dx      = cr3bpEomControl(u,p.crp,t,at)
    dλ      = cr3bpCostateEom(u,p, γ)
    usqr    = γ^2
    du = @SVector [dx[1], dx[2], dx[3], dx[4], dx[5], dx[6], -γ*tMaxSc / cSc,
                   dλ[1], dλ[2], dλ[3], dλ[4], dλ[5], dλ[6], dλ[7], usqr]

    return du
    end
end

function cr3bpEomIndirectEnergyIntegral!(du::AbstractArray, u::AbstractArray, p::AbstractCR3BPIndirectParams, t)
    @inbounds begin
    # Get requirements 
    TU  = p.crp.TU
    LU  = p.crp.LU
    VU  = p.crp.VU 
    MU  = p.crp.MU 
    isp = p.sp.isp 
    c   = p.sp.c

    # Scale requirements
    tMaxSc = p.sp.tMax * TU^2 / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]^2 + u[12]^2 + u[13]^2)
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p.utype, p.ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Compute dynamics
    cr3bpEomControl!(view(du,1:6), u, p.crp, t, at)
    du[7]   = -γ*tMaxSc / cSc
    GVec    = cr3bpCostateEom!(view(du,8:14), u, p, γ)
    du[15]  = γ^2
    end
    return GVec
end

# Utility Functions 
function computeS(x::AbstractVector, λv, cSc)
    @inbounds begin
        return -λv*cSc / x[7] - x[14] + 1.0
    end
end

function computeU(S, utype, ϵ)
    if utype == 2
        return 1.0
    elseif utype == 0
        return 0.0
    else
        return (ϵ - S) / (2.0*ϵ)
    end
end





































