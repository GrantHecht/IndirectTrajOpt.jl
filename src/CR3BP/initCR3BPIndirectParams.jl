
# Instantiates CR3BPIndirectParams struct for a given scenario
function initCR3BPIndirectParams(scenario::String)
    if scenario == "Low Thrust 10 CR3BP"
        sp = SimpleSpacecraft(1500.0, 1000.0, 10.0, 3000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          1500.0, 384400)
        ps = CR3BPIndirectParams(sp, cr3, 1.0)
    elseif scenario ==  "Low Thrust 1 CR3BP"
        sp = SimpleSpacecraft(1500.0, 1000.0, 1.0, 3000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          1500.0, 384400)
        ps = CR3BPIndirectParams(sp, cr3, 1.0)
    elseif scenario == "Low Thrust 1.5 CR3BP"
        sp = SimpleSpacecraft(2000.0, 1000.0, 1.5, 2000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          2000.0, 384400)
        ps = CR3BPIndirectParams(sp, cr3, 1.0) 
    elseif scenario == "Dawn CR3BP"
    	sp = SimpleSpacecraft(1218.0, 1218.0, 92.0e-3, 3200.0)
	cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          2000.0, 384400)
        ps = CR3BPIndirectParams(sp, cr3, 1.0) 
    else
        error("Scenario is not implemented!")
    end

    return ps
end

# Instantiates CR3BPIndirectWithSTM struct for a given scenario
function initCR3BPIndirectWithSTMParams(scenario::String)
    if scenario == "Low Thrust 10 CR3BP"
        sp = SimpleSpacecraft(1500.0, 1000.0, 10.0, 3000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          1500.0, 384400)
        ps = CR3BPIndirectWithSTMParams(sp, cr3, 1.0)
    elseif scenario == "Low Thrust 1 CR3BP"
        sp = SimpleSpacecraft(1500.0, 1000.0, 1.0, 3000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          1500.0, 384400)
        ps = CR3BPIndirectWithSTMParams(sp, cr3, 1.0)
    elseif scenario == "Low Thrust 1.5 CR3BP"
        sp = SimpleSpacecraft(2000.0, 1000.0, 1.5, 2000.0)
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          2000.0, 384400)
        ps = CR3BPIndirectWithSTMParams(sp, cr3, 1.0) 
    elseif scenario == "Dawn CR3BP"
    	sp = SimpleSpacecraft(1218.0, 1218.0, 92.0e-3, 3200.0)
	cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          2000.0, 384400)
        ps = CR3BPIndirectWithSTMParams(sp, cr3, 1.0) 
    else
        error("Scenario is not implemented!")
    end

    return ps
end
