
function cr3bpEomControlMTKSys(μ::AbstractFloat)
    @parameters atx aty atz
    @variables t x(t) y(t) z(t) vx(t) vy(t) vz(t)
    D = Differential(t)

    # Compute r1 and r2
    r1 = sqrt((x+μ)^2 + y^2 + z^2)
    r2 = sqrt((x+μ-1)^2 + y^2 + z^2)

    # Compute g 
    gx = x - (1 - μ)*(x + μ)/r1^3 - μ*(x + μ - 1)/r2^3
    gy = y - (1 - μ)*y/r1^3 - μ*y/r2^3
    gz = -(1 - μ)*z/r1^3 - μ*z/r2^3
    
    # Compute h
    hx = 2*vy
    hy = -2*vx
    
    eqs = [D(x) ~ vx, 
           D(y) ~ vy, 
           D(z) ~ vz, 
           D(vx) ~ gx + hx + atx,
           D(vy) ~ gy + hy + aty,
           D(vz) ~ gz + atz]

    return ODESystem(eqs)
end

function cr3bpEomNoControlMTKSys(μ::AbstractFloat)
    _cr3bpEomWControl = cr3bpEomControlMTKSys(μ)
    (atx, aty, atz) = parameters(_cr3bpEomWControl)
    eqs = [eq.lhs ~ simplify(
        substitute(substitute(substitute(
            eq.rhs, atx => 0), aty => 0), atz => 0)) for 
            eq in equations(_cr3bpEomWControl)]
    return ODESystem(eqs)
end

function cr3bpEOMIndirect(μ::AbstractFloat,    TU::AbstractFloat,
                          LU::AbstractFloat,   MU::AbstractFloat,
                          tMax::AbstractFloat, c::AbstractFloat)
    @parameters u 
    @variables t x(t) y(t) z(t) vx(t) vy(t) vz(t) m(t) 
    @variables λx(t) λy(t) λz(t) λvx(t) λvy(t) λvz(t) λm(t) H(t)
    D    = Differential(t)

    # Compute r1 and r2
    r1 = sqrt((x+μ)^2 + y^2 + z^2)
    r2 = sqrt((x+μ-1)^2 + y^2 + z^2)

    # Compute g 
    gx = x - (1 - μ)*(x + μ)/r1^3 - μ*(x + μ - 1)/r2^3
    gy = y - (1 - μ)*y/r1^3 - μ*y/r2^3
    gz = -(1 - μ)*z/r1^3 - μ*z/r2^3
    
    # Compute h
    hx = 2*vy
    hy = -2*vx

    # Compute G 
    G11     = 1 - (1 - μ)/r1^3 + 3*(1 - μ)*(x + μ)^2/r1^5 - 
                μ/r2^3 + 3*μ*(x + μ - 1)^2/r2^5;
    G22     = 1 - (1 - μ)/r1^3 + 3*(1 - μ)*y^2/r1^5 - 
                μ/r2^3 + 3*μ*y^2/r2^5;
    G33     = -(1 - μ)/r1^3 + 3*(1 - μ)*z^2/r1^5 - 
                μ/r2^3 + 3*μ*z^2/r2^5;
    G12     = 3*(1 - μ)*(x + μ)*y/r1^5 + 
                3*μ*(x + μ - 1)*y/r2^5;
    G13     = 3*(1 - μ)*(x + μ)*z/r1^5 + 
                3*μ*(x + μ - 1)*z/r2^5;
    G23     = 3*(1 - μ)*y*z/r1^5 + 
                3*μ*y*z/r2^5;

    # Compute scaled max thrust and exhaust velocity
    tMaxSc = tMax * TU * TU / (MU*LU*1000.0)
    cSc = c * TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(λvx^2 + λvy^2 + λvz^2)
    αx = -λvx / λv 
    αy = -λvy / λv 
    αz = -λvz / λv

    # Compute thrust acceleration
    atx = u*tMaxSc*αx / m
    aty = u*tMaxSc*αy / m
    atz = u*tMaxSc*αz / m

    # Dynamics
    dx   = vx
    dy   = vy
    dz   = vz
    ddx  = gx + hx + atx
    ddy  = gy + hy + aty
    ddz  = gz + atz
    dm   = -u*tMaxSc / cSc
    dλx  = -G11*λvx - G12*λvy - G13*λvz
    dλy  = -G12*λvx - G22*λvy - G23*λvz
    dλz  = -G13*λvx - G23*λvy - G33*λvz
    dλvx = -λx + 2.0*λvy
    dλvy = -λy - 2.0*λvx
    dλvz = -λz
    dλm  = -λv*u*tMaxSc / (m^2)

    # Equations
    eqs = [D(x)   ~ dx,
           D(y)   ~ dy,
           D(z)   ~ dz,
           D(vx)  ~ ddx,
           D(vy)  ~ ddy,
           D(vz)  ~ ddz,
           D(m)   ~ dm,
           D(λx)  ~ dλx,
           D(λy)  ~ dλy,
           D(λz)  ~ dλz,
           D(λvx) ~ dλvx,
           D(λvy) ~ dλvy,
           D(λvz) ~ dλvz,
           D(λm)  ~ dλm]
    
    return ODESystem(eqs, t)
end