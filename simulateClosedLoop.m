function data = simulateClosedLoop(x0, theta, paramsPlant, params, mpcOpt, simOpt)

    if ~isfield(simOpt, 'zTarget'), simOpt.zTarget = 0.8; end

    maxSteps = simOpt.M;
    x = x0;
    nx = 3;
    xTrj = zeros(nx, maxSteps+1); xTrj(:,1) = x0;
    uTrj = zeros(1, maxSteps);
    VTTrj = zeros(1,maxSteps);

    violated = struct('VT', false, 'T', false);
    t = (0:maxSteps)*paramsPlant.dt;

    for k = 1:maxSteps
        [uSeq, ~, ~] = mpcController(x, theta, params, mpcOpt);
        u = uSeq(1);
        [xNext, VT] = ecmDynamics(x, u, paramsPlant);

        xTrj(:,k+1) = xNext;
        uTrj(k)     = u;
        VTTrj(k)    = VT;

        if VT > mpcOpt.VT_max || VT < mpcOpt.VT_min
            violated.VT = true;
        end
        if xNext(3) > mpcOpt.T_max
            violated.T = true;
        end

        x = xNext;
        if x(1) >= simOpt.zTarget
            break;
        end
    end

    k_end = k;
    data.xTrj = xTrj(:,1:k_end+1);
    data.uTrj = uTrj(1:k_end);
    data.VT    = VTTrj(1:k_end);
    data.t     = t(1:k_end+1);
    data.violations = violated;
    data.tCharge   = data.t(end);
end