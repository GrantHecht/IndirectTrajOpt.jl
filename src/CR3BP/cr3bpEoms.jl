
# CR3BP Parameters
struct CR3BPParams 

    # Primary and Secondary Body Mass
    m1::Float64
    m2::Float64

    # Radius of Primary and Secondary Bodies
    R1::Float64
    R2::Float64

    # Total Mass
    mtot::Float64

    # Gravitational Parameters
    μ::Float64

    # Length Unit
    LU::Float64

    # Time Unit
    TU::Float64

    # Speed Unit
    VU::Float64

    # Mass Unit
    MU::Float64

    # Gravitational Constant
    G::Float64

    function CR3BPParams(m1,m2,R1,R2,MU,r12) 
        G       = 6.673e-20
        mtot    = m1 + m2
        μ       = m2 / (m1 + m2)
        LU      = r12 
        TU      = 1.0 / sqrt(G * (m1 + m2) / (LU^3))
        VU      = LU / TU 
        MU      = MU 
        new(m1,m2,R1,R2,mtot,μ,LU,TU,VU,MU,G)
    end
end

# CR3BP Eq. of Motion w/o control
function cr3bpEomNoControl(u::AbstractArray, p::CR3BPParams, t)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    return @SVector [u[4], u[5], u[6],
                     gx + hx,
                     gy + hy,
                     gz]
end
function cr3bpEomNoControl!(du::AbstractArray,u::AbstractArray,p::CR3BPParams,t)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Set du 
    @. du[1:3] = view(u,4:6)
    du[4] = gx + hx 
    du[5] = gy + hy
    du[6] = gz
end

# CR3BP Eq. of Motion w/ control
function cr3bpEomControl(u::AbstractArray,p::CR3BPParams,t,
                         a::AbstractArray)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = inv(r1^3)
    invR23 = inv(r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Compute and return dynamics
    dx  = u[4]
    dy  = u[5]
    dz  = u[6]
    ddx = gx + hx + a[1]
    ddy = gy + hy + a[2]
    ddz = gz + a[3]

    return @SVector [dx,dy,dz,ddx,ddy,ddz]
end

function cr3bpEomControl!(du::AbstractArray,u::AbstractArray,p::CR3BPParams,t,
                          a::AbstractArray)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Set du 
    @inbounds for i in 1:3
        du[i] = u[i + 3]
    end
    du[4] = gx + hx + a[1]
    du[5] = gy + hy + a[2]
    du[6] = gz + a[3]
end
  