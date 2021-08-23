
# Instantiates CR3BPIndirectParams struct for a given scenario
function initCR3BPIndirectParams(scenario::String)
    if scenario == "Low Thrust CR3BP"

        # Initialize Spacecraft
        sp = SimpleSpacecraft(1500.0, 1000.0, 10.0, 3000.0)

        # Initialize CR3BP params 
        cr3 = CR3BPParams(5.9742e24, 7.3483e22,
                          6378.137, 1738.0,
                          1500.0, 384400)

        # Initialize Indirect CR3BP Optimization struct 
        ps = CR3BPIndirectParams(sp, cr3, 1.0)
    else
        error("Scenario is not implemented!")
    end
end