function thetaNext = proposeNextTheta(gpJ, gpG_list, thetaLog, JLog, GiLog, opt)
    Ncand = opt.Ncand;
    lb    = opt.lb;
    ub    = opt.ub;
    kappa = opt.kappa;   % for objective
    beta  = opt.beta;    % for constraints

    tau   = 3.0;         % barrier weight

    nTheta = numel(lb);

    % Random candidate thetas in [lb, ub]
    ThetaCand = bsxfun(@plus, lb, bsxfun(@times, (ub-lb), rand(Ncand, nTheta)));

    acqVals = -inf(Ncand,1);

    for j = 1:Ncand
        th = ThetaCand(j,:);

        [muJ, sJ] = predict(gpJ, th);

        % Acquisition for objective (we minimize J -> LCB in J -> negative)
        acq0 = -(muJ - kappa * sJ);

        acqBarrier = 0;
        for i = 1:numel(gpG_list)
            [muGi, sGi] = predict(gpG_list{i}, th);
            lcbGi = muGi - beta * sGi;

            if lcbGi <= 0
                acqBarrier = -inf;
                break;
            else
                acqBarrier = acqBarrier + log(lcbGi);
            end
        end

        acqVals(j) = acq0 + tau * acqBarrier;
    end

    [~, idx] = max(acqVals);
    thetaNext = ThetaCand(idx,:).';
end
