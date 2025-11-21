function metrics = runRLEpisode(x0, Q, SOCBins, TBins, ...
                                  paramsEnv, VTMax, VTMin, ...
                                  TMax, IMax, zTarget, ...
                                  maxSteps, actionValues)

    dt = paramsEnv.dt;


    [~, VT0] = batteryStep(x0, 0, paramsEnv);
    s = [x0(1); VT0; x0(2); x0(3)];

    zHist  = zeros(1, maxSteps+1);
    THist  = zeros(1, maxSteps+1);
    VTHist = zeros(1, maxSteps+1);
    IHist  = zeros(1, maxSteps);
    tVec   = (0:maxSteps)*dt;


    zHist(1)  = x0(1);
    THist(1)  = x0(3);
    VTHist(1) = VT0;


    chargingTime = inf;
    violated = false;


    for k = 1:maxSteps
        s_idx       = discretizeState(s, SOCBins, TBins);
        [~, a_idx]  = max(Q(s_idx,:));
        I           = actionValues(a_idx);


        % Step environment
        x = [s(1); s(3); s(4)];
        [xNext, VTNext] = batteryStep(x, I, paramsEnv);
        zNext  = xNext(1);
        v1Next = xNext(2);
        TNext  = xNext(3);


        IHist(k)      = I;
        zHist(k+1)    = zNext;
        THist(k+1)    = TNext;
        VTHist(k+1)   = VTNext;

        s = [zNext; VTNext; v1Next; TNext];


        % Safety check
        if (VTNext > VTMax) || (VTNext < VTMin) || ...
           (TNext > TMax)   || (abs(I) > IMax)
            violated     = true;
            chargingTime = tVec(k);
            kEnd = k;
            break;
        end


        % Target check
        if zNext >= zTarget && isinf(chargingTime)
            chargingTime = tVec(k);
            kEnd = k;
            break;
        end

        kEnd = k;
    end

    zHist  = zHist(1:kEnd+1);
    THist  = THist(1:kEnd+1);
    VTHist = VTHist(1:kEnd+1);
    IHist  = IHist(1:kEnd);
    t       = tVec(1:kEnd+1);

    metrics = struct();
    metrics.chargingTime = chargingTime;
    metrics.violated = violated;
    metrics.zHist = zHist;
    metrics.THist = THist;
    metrics.VT_Hist = VTHist;
    metrics.IHist = IHist;
    metrics.t = t;
end