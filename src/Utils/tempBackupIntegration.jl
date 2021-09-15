
function tempIntegrate!(x, tspan, ps; avoidErr = false)
    h0 = 0.5 / ps.crp.TU
    maxh = 10.0 / ps.crp.TU

    # Set utype 
    cSc = ps.sp.isp*9.81*ps.crp.TU / (ps.crp.LU*1000.0)
    λv = norm(view(x,11:13))
    S = computeS(x, λv, cSc)
    if S > ps.ϵ
        ps.utype = 0
    elseif S < -ps.ϵ
        ps.utype = 2
    else
        ps.utype = 1
    end

    # Initialize integrator
    STM = true
    if length(x) == 210
        integ = ATO_Integraion.initODE((dz,z,p,t)->cr3bpEomIndirectWithSTM!(dz,z,p,t),
                                       x,
                                       tspan,
                                       ATO_Integraion.Vern9ME();
                                       ps = ps,
                                       RelTol = 1e-14,
                                       AbsTol = 1e-14,
                                       maxh = maxh)
        x0sm = Vector{Float64}(undef, 14)
        x0sm .= view(x,1:14)
        integ_sw = ATO_Integraion.initODE((dz,z,p,t)->cr3bpEomIndirect!(dz,z,p,t),
                                          x0sm,
                                          tspan,
                                          ATO_Integraion.Vern9ME();
                                          ps = ps,
                                          RelTol = 1e-14,
                                          AbsTol = 1e-14,
                                          maxh = maxh)
    elseif length(x) == 14
        STM = false
        integ = ATO_Integraion.initOde((dz,z,p,t)->cr3bpEomIndirect!(dz,z,p,t),
                                       x,
                                       tspan,
                                       ATO_Integraion.Vern9ME();
                                       ps = ps,
                                       RelTol = 1e-14,
                                       AbsTol = 1e-14,
                                       maxh = maxh)
    end

    # Begin integration loop
    done = false
    switching = false
    reducing = false
    sw_cnt = 0
    while !done

        # Incriment Switching Counter
        sw_cnt += 1

        # Take Step
        ATO_Integration.step!(integ)

        # Check for Switching
        switching, reducing = Check_4_Switching(integ.ukp1, integ.ps)

        # Detect Switching Time
        if switching && sw_cnt > ss_iters # Switching
            if STM
                Update_IntegSW!(integ_sw, integ)
                tsw = Detect_t_Switch!(integ_sw, mn_iters, mr_iters, s_tol, avoid_err)
            else
                tsw = Detect_t_Switch!(integ, mn_iters, mr_iters, s_tol, avoid_err)
            end
            if tsw < 0 || tsw == integ.tk
                switching = false
                reducing = false

                # Make sure utype is corect, reduce step size, and restep
                _Correct_Utype!(integ)
                integ.h = h0
                ATO_Integration.step!(integ)
            end
        end

        # Handle Switching or Reducing
        if switching && sw_cnt > ss_iters # Switching

            # Reset Switching Counter
            sw_cnt = 0

            # Detect Switching Time
            #if STM
            #    Update_IntegSW!(integ_sw, integ)
            #    tsw = Detect_t_Switch!(integ_sw, mn_iters, mr_iters, s_tol, avoid_err)
            #else
            #    tsw = Detect_t_Switch!(integ, mn_iters, mr_iters, s_tol, avoid_err)
            #end

            # Step to tsw
            ATO_Integration.step!(integ, tsw - integ.tk)

            # Update Utype
            _Update_Utype!(integ)

            # Accept Step
            ATO_Integration.accept_step!(integ)

            # Update Step Size
            integ.h = integ.h * reduce_fact
            ATO_Integration.init_h!(integ)

            # Propagate STM over Discontinuity
            if params.ϵ == 0.0 && STM
                _Prop_STM_Discont!(integ)
            end

        elseif reducing # Reducing
            integ.h = integ.h * reduce_fact

        else # Continue to Next Step W/O Doing Anything
            ATO_Integration.accept_step!(integ)
        end

        # Handle Closeness to tf
        if integ.tk >= integ.tf
            done = true
        elseif integ.tk + integ.h >= integ.tf
            integ.h = integ.tf - integ.tk
        end
    end

    # Update x before returning
    @inbounds for i in 1:length(x)
        x[i] = integ.uk[i]
    end
end

function Check_4_Switching(x::AbstractArray, params)

    # Set Bools
    switching = false
    reducing = false

    # Compute Skp1
    Skp1 = compute_S(x, params)

    # Check for Switching
    if params.utype == 2 # Thrust is On
        if Skp1 < -params.ϵ # No Switching
            switching = false
        elseif Skp1 > -params.ϵ && Skp1 < params.ϵ # Switching to Medium
            switching = true
        else # Reducing Step Size or Switching
            if params.ϵ == 0.0
                switching = true
            else
                reducing = true
            end
        end
    elseif params.utype == 0 # Thrust is Off

        if Skp1 > params.ϵ # No Switching
            switching = false
        elseif Skp1 > -params.ϵ && Skp1 < params.ϵ # Switching to Medium
            switching = true
        else # Reducing Step Size or Switching
            if params.ϵ == 0.0
                switching = true
            else
                reducing = true
            end
        end
    else # Medium Thrust
        if !(Skp1 > -params.ϵ && Skp1 < params.ϵ) # Switching On or Off
            switching = true
        end
    end

    return switching, reducing
end

function Update_IntegSW!(integ_sw::T,integ::U) where {T,U<:ATO_Integration.ODEInt}
    integ_sw.tk = integ.tk
    integ_sw.tkp1 = integ.tkp1
    integ_sw.h = integ.h
    integ_sw.qold = integ_sw.qold
    integ_sw.ps.utype = integ.ps.utype
    @inbounds for i in 1:14
        integ_sw.uk[i] = integ.uk[i]
        integ_sw.ukp1[i] = integ.ukp1[i]
    end
end

function Detect_t_Switch!(integ::ATO_Integration.ODEInt, mn_iters, mr_iters, s_tol, avoid_err)

    # Initialize tsw
    tkp1 = integ.tkp1
    tsw = integ.tk

    # Initialize Switching Function
    Sjp1 = 0.0
    if integ.ps.utype == 1
        Sjp1 = compute_S(integ.ukp1, integ.ps)
    end
    _fsw(tsw) = _Switch_Function!(tsw, Sjp1, integ)

    # Begin Newton's Itterations
    stop = false
    succ = false
    niters = 0
    while !stop
        # Incriment Iteration Counter
        niters += 1

        # Compute Switching Function
        fsw = _fsw(tsw)

        # Compute 1/Sdot(tsw)
        if tsw > integ.tk
            dSinv = _Compute_dSinv(integ.ukp1, integ.ps)
        else
            dSinv = _Compute_dSinv(integ.uk, integ.ps)
        end

        # Update tsw
        tsw -= fsw*dSinv

        # Check for proper convergence
        if tsw < integ.tk || tsw > integ.tkp1 || niters > mn_iters
            stop = true
        elseif abs(fsw) < s_tol
            stop = true
            succ = true
        end
    end

    # Begin Regula-Falsi is Newtons Failed
    if !succ
        # Initialize
        stop = false
        niters = 0
        a = integ.tk
        b = tkp1
        fa = _fsw(a)
        fb = _fsw(b)

        # Check for root bracketing
        if fa * fb > 0
            if avoid_err
                #if abs(fa) < abs(fb)
                #    tsw = a
                #else
                #    tsw = b
                #end
                tsw = -1.0
                stop = true
            else
                throw(ErrorException("Root not bracketed for Regula-Falsi switching detection!"))
            end
        end

        # Regula Falsi Loop
        while !stop

            # Incriment Iteration Counter
            niters += 1

            # Compute c
            c = b - (fb*(b - a))/(fb - fa)
            fc = _fsw(c)

            # Update
            if fa*fc < 0
                b = c
                fb = fc
            else
                a = c
                fa = fc
            end

            # Update tsw and check for proper convergence
            if abs(fa) < abs(fb)
                fsw = fa
                tsw = a
            else
                fsw = fb
                tsw = b
            end
            if abs(fsw) < s_tol
                stop = true
                succ = true
            elseif niters > mr_iters
                stop = true
            end
        end
    end

    return tsw
end

function _Switch_Function!(tsw, Sjp1, integ::ATO_Integration.ODEInt)

    # Determine Switching Direction
    ϵ = -integ.ps.ϵ
    if integ.ps.utype == 2 || (integ.ps.utype == 1 && Sjp1 < -integ.ps.ϵ)
        ϵ *= -1.0
    end

    # Take Step
    h = tsw - integ.tk
    if h > 0
        ATO_Integration.step!(integ, h)
        Stsw = compute_S(integ.ukp1, integ.ps)
    else
        Stsw = compute_S(integ.uk, integ.ps)
    end

    # Compute Switching Function
    return Stsw + ϵ
end

function _Compute_dSinv(x::AbstractVector, params)
    λ_v = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    c_nsc = params.sp.Isp * 9.81 * params.eom.TU / (params.eom.LU * 1000)

    # Compute 1/Sdot(tsw)
    dSinv = (x[8] - 2*x[12])*x[11] + (x[9] + 2*x[11])*x[12] + x[10]*x[13]
    dSinv = λ_v*x[7]/(c_nsc*dSinv)

    return dSinv
end

function _Update_Utype!(integ::ATO_Integration.ODEInt)
    if (integ.ps.utype == 2 || integ.ps.utype == 0) && integ.ps.ϵ != 0.0
        integ.ps.utype = 1
    elseif integ.ps.ϵ != 0.0
        Stsw = compute_S(integ.ukp1, integ.ps)
        if abs(Stsw + integ.ps.ϵ) < abs(Stsw - integ.ps.ϵ)
            integ.ps.utype = 2
        else
            integ.ps.utype = 0
        end
    else
        if integ.ps.utype == 2
            integ.ps.utype = 0
        else
            integ.ps.utype = 2
        end
    end
end

function _Correct_Utype!(integ::ATO_Integration.ODEInt)
    S = compute_S(integ.uk, integ.ps)
    if S < -integ.ps.ϵ
        integ.ps.utype = 2
    elseif S > integ.ps.ϵ
        integ.ps.utype = 0
    else
        integ.ps.utype = 1
    end
end

function _Prop_STM_Discont!(integ::ATO_Integration.ODEInt)

    # Revert Utype of Parameters Copy to Previous Value
    integ.ps.utype == 2 ? integ.ps.utype = 0 : integ.ps.utype = 2

    # Initialize Requirements
    temp = MVector{14, Float64}(undef)
    dy⁻ = MVector{14,Float64}(undef)
    dy⁺ = MVector{14,Float64}(undef)

    # Compute dy⁻
    @inbounds for i in 1:14; temp[i] = integ.uk[i]; end
    CR3BPIMF_EOM!(dy⁻, temp, integ.ps, 0.0)

    # Compute dy⁺
    integ.ps.utype == 2 ? integ.ps.utype = 0 : integ.ps.utype = 2
    @inbounds for i in 1:14; temp[i] = integ.uk[i]; end
    CR3BPIMF_EOM!(dy⁺, temp, integ.ps, 0.0)

    # Compute ∂S and dSinv
    λ_v = sqrt(integ.uk[11]^2 + integ.uk[12]^2 + integ.uk[13]^2)
    c_nsc = integ.ps.sp.Isp * 9.81 * integ.ps.eom.TU / (integ.ps.eom.LU * 1000)

    @views @. temp .= [0.0;0.0;0.0;0.0;0.0;0.0;
                       λ_v*c_nsc / (integ.uk[7]^2);
                       0.0;0.0;0.0;
                       integ.uk[11:13] * (-c_nsc / (λ_v*integ.uk[7]));
                       -1.0]
    dSinv = (integ.uk[8] - 2*integ.uk[12])*integ.uk[11] +
            (integ.uk[9] + 2*integ.uk[11])*integ.uk[12] +
            (integ.uk[10])*integ.uk[13]
    dSinv = λ_v*integ.uk[7] / (c_nsc*dSinv)

    # Compute Φ and Ψ
    uSTM = reshape(view(integ.uk, 15:210), 14, 14)
    Φ = SizedMatrix{14,14,Float64}(undef)
    Ψ = SizedMatrix{14,14,Float64}(undef)
    Threads.@threads for row in 1:14
        for col in 1:14
            Φ[row, col] = uSTM[row, col]
            Ψ[row, col] = I[row,col] + (dy⁺[row] - dy⁻[row]) * temp[col]*dSinv
        end
    end

    # Multiplying Ψ and Φ
    Threads.@threads for row in 1:14
        for col in 1:14
            sum = 0.0
            for k in 1:14
                sum += Ψ[row, k] * Φ[k, col]
            end
            uSTM[row, col] = sum
        end
    end
end

function compute_S(x::AbstractVector, params)
    λ_v = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    c_nsc   = params.sp.Isp * 9.81 * params.eom.TU / (params.eom.LU * 1000)

    return (-λ_v * c_nsc / x[7] - x[14]  + 1.0)
end