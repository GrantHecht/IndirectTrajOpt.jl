
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