function [sNext, r, done] = rlStep(s, a_idx, paramsEnv, ...
                                    VTMax, VTMin, TMax, IMax, ...
                                    zTarget, actionValues)

    z  = s(1);
    VT = s(2);
    v1 = s(3);
    T  = s(4);

    I = actionValues(a_idx);

    x = [z; v1; T];
    [xNext, VTNext] = batteryStep(x, I, paramsEnv);

    zNext  = xNext(1);
    v1Next = xNext(2);
    TNext  = xNext(3);

    sNext = [zNext; VTNext; v1Next; TNext];


    dSOC = zNext - z;

    r = 100*dSOC - 0.1;


    violated = (VTNext > VTMax) || (VTNext < VTMin) || (TNext  > TMax)  || (abs(I) > IMax);

    if violated
        r    = r - 1000;
        done = true;
        return;
    end

    if zNext >= zTarget
        r    = r + 500;
        done = true;
        return;
    end

    done = false;
end
