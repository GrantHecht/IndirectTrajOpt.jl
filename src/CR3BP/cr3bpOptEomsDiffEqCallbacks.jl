
function cr3bpEomsCondition(out, x, t, integrator)

    # Get Requirements
    TU  = integrator.p.crp.TU
    LU  = integrator.p.crp.LU
    MU  = integrator.p.crp.MU 
    R1  = integrator.p.crp.R1
    R2  = integrator.p.crp.R2
    μ   = integrator.p.crp.μ
    isp = integrator.p.sp.isp
    m0  = integrator.p.sp.initMass
    mp  = integrator.p.sp.initProp

    ϵ   = integrator.p.ϵ
    utype = integrator.p.utype

    # Compute Scaling
    cSc = isp*9.81*TU / (LU * 1000.0)
    
    # Compute required
    λv = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    S  = computeS(x, λv, cSc)
    xpμ = x[1] + μ
    r1  = sqrt(xpμ^2 + x[2]^2 + x[3]^2)
    r2  = sqrt((xpμ - 1)^2 + x[2]^2 + x[3]^2) 

    # Switching condition 
    if ϵ != 0.0
        if  utype == 0
            out[1] = S - ϵ
        elseif utype == 2
            out[1] = S + ϵ
        else
            out[1] = abs(S) - ϵ
        end
    else
        out[1] = S
    end

    # Termination conditions 
    out[2] = x[7]*MU - (m0 - mp)
    out[3] = r1*LU - R1 
    out[4] = r2*LU - R2 
end

function cr3bpEomsAffect!(integrator, idx)

    # Get Requirements
    TU  = integrator.p.crp.TU
    LU  = integrator.p.crp.LU
    isp = integrator.p.sp.isp

    ϵ   = integrator.p.ϵ
    utype = integrator.p.utype

    # Switching affect 
    if idx == 1
        if ϵ != 0.0 
            if utype != 1
                integrator.p.utype = 1
            else
                cSc = isp*9.81*TU / (LU*1000.0)
                λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
                S = computeS(integrator.u, λv, cSc)
                if S < 0.0
                    integrator.p.utype = 2
                else
                    integrator.p.utype = 0
                end
            end

        else # Need to prop STM across switching 
            if utype == 0
                integrator.p.utype = 2
            else
                integrator.p.utype = 0
            end
        end

        # Reset step 
        auto_dt_reset!(integrator)

    # Termination affect
    else
        println("Terminating!")
        terminate!(integrator)
    end
end

function cr3bpEomsConditionNoTerm(x, t, integrator)

    # Get Requirements
    TU  = integrator.p.crp.TU
    LU  = integrator.p.crp.LU
    isp = integrator.p.sp.isp

    ϵ   = integrator.p.ϵ
    utype = integrator.p.utype

    # Compute Scaling
    cSc = isp*9.81*TU / (LU * 1000.0)
    
    # Compute required
    λv = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    S  = computeS(x, λv, cSc)

    # Switching condition 
    out = 0.0
    if ϵ != 0.0
        if  utype == 0
            out = S - ϵ
        elseif utype == 2
            out = S + ϵ
        else
            out = abs(S) - ϵ
        end
    else
        out = S
    end

    return out
end

function cr3bpEomsAffectNoTerm!(integrator)

    # Get Requirements
    TU  = integrator.p.crp.TU
    LU  = integrator.p.crp.LU
    isp = integrator.p.sp.isp

    ϵ   = integrator.p.ϵ
    utype = integrator.p.utype

    # Switching affect 
    if ϵ != 0.0 
        if utype != 1
            integrator.p.utype = 1
        else
            cSc = isp*9.81*TU / (LU*1000.0)
            λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
            S = computeS(integrator.u, λv, cSc)
            if S < 0.0
                integrator.p.utype = 2
            else
                integrator.p.utype = 0
            end
        end

    else # Need to prop STM across switching 
        if utype == 0
            integrator.p.utype = 2
        else
            integrator.p.utype = 0
        end
    end

    # Reset step 
    auto_dt_reset!(integrator)
end

function cr3bpEomsAffectNoTermWithSTM!(integrator)

    # Get Requirements
    TU  = integrator.p.crp.TU
    LU  = integrator.p.crp.LU
    isp = integrator.p.sp.isp

    ϵ   = integrator.p.ϵ
    utype = integrator.p.utype

    # Switching affect 
    if ϵ != 0.0 
        if utype != 1
            integrator.p.utype = 1
        else
            cSc = isp*9.81*TU / (LU*1000.0)
            λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
            S = computeS(integrator.u, λv, cSc)
            if S < 0.0
                integrator.p.utype = 2
            else
                integrator.p.utype = 0
            end
        end
    else 
        propSTM!(integrator.u, integrator.p)
    end

    # Reset step 
    auto_dt_reset!(integrator)
end

function propSTM!(u::AbstractVector, ps::CR3BPIndirectWithSTMParams)
    @inbounds begin
        # Get requirements 
        TU  = p.crp.TU
        LU  = p.crp.LU
        MU  = p.crp.MU 
        μ   = p.crp.μ
        c   = p.sp.c

        # Scale Requirements
        tMaxSc = p.sp.tMax * TU * TU / (MU*LU*1000.0)
        cSc = c*TU / (LU*1000.0)

        # Compute Requirements 
        λv = sqrt(u[11]^2 + u[12]^2 + u[13]^2)

        # Get state/co-state derivatives before updating utype
        dy⁻ = cr3bpEomIndirect(view(u, 1:14), ps, 0.0)

        # Update utype 
        ps.utype == 0 ? ps.utype = 2 : ps.utype = 0

        # Get state/co-state derivatives after updating utype
        dy⁺ = cr3bpEomIndirect(view(u, 1:14), ps, 0.0)

        # Compute dy difference
        dyDiff = @SVector [dy⁺[1]  - dy⁻[1],
                           dy⁺[2]  - dy⁻[2],
                           dy⁺[3]  - dy⁻[3],
                           dy⁺[4]  - dy⁻[4],
                           dy⁺[5]  - dy⁻[5],
                           dy⁺[6]  - dy⁻[6],
                           dy⁺[7]  - dy⁻[7],
                           dy⁺[8]  - dy⁻[8],
                           dy⁺[9]  - dy⁻[9],
                           dy⁺[10] - dy⁻[10],
                           dy⁺[11] - dy⁻[11],
                           dy⁺[12] - dy⁻[12],
                           dy⁺[13] - dy⁻[13],
                           dy⁺[14] - dy⁻[14]]

        # Compute 1 / Sdot and ∂S/∂y / Sdot
        dSInv   = (u[8] - 2.0*u[12])*u[11] + (u[9] + 2.0*u[11])*u[12] + u[10]*u[13]
        temp    = -cSc*dSInv / (u[7]*λv)
        ∂SdSInv = @SVector [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                            λv*cSc*dSInv / (u[7]^2), 
                            0.0, 0.0, 0.0,
                            u[11]*temp, 
                            u[12]*temp, 
                            u[13]*temp, 
                            -1.0*dSInv]

        # Compute Ψ using allocated jac matrix in parameters 
        Ψ = ps.m1; Ψ .= 0.0
        for i in 1:14
            Ψ[i,i] = 1.0
        end
        mul!(Ψ, dyDiff, transpose(∂SdSInv), 1.0, 1.0)

        # Propagate STM 
        Φ⁺ = reshape(view(u, 15:210), (14,14))
        Φ⁻ = ps.m2; Φ⁻ .= Φ⁺
        mul!(Φ⁺, Ψ, Φ⁻)
    end
    return nothing
end