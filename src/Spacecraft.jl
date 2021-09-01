# Spacecraft Models
abstract type Spacecraft end

# Simple Spacecraft Model
struct SimpleSpacecraft <: Spacecraft

    # Initial mass 
    initMass::Float64

    # Initial propellant mass 
    initProp::Float64 

    # Max thrust 
    tMax::Float64

    # Specific impulse 
    isp::Float64

    # Exaust velocity
    c::Float64

    function SimpleSpacecraft(initMass, initProp, tMax, isp)
        c = 9.81*isp 
        new(initMass, initProp, tMax, isp, c)
    end
end