function [xNext, VT] = batteryStep(x, I, params)
    [xNext, VT] = ecmDynamics(x, I, params);
end
