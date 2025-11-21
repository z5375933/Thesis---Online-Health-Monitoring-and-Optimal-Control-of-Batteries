function [J, VTPeak, TPeak, IPeak, data, metrics] = evaluateTheta(theta, x0, paramsPlant, params, mpcOpt, simOpt)
   data = simulateClosedLoop(x0, theta, paramsPlant, params, mpcOpt, simOpt);

   SOCTraj = data.xTrj(1,:);
   TTraj   = data.xTrj(3,:);
   VTraj  = data.VT;
   uTraj   = data.uTrj;

   SOCTarget = simOpt.z_target;
   idx_reached = find(SOCTraj >= SOCTarget, 1, 'first');
   if isempty(idx_reached)
       chargingTime = inf;
   else
       chargingTime = data.t(idx_reached);
   end

   dt = paramsPlant.dt;
   if isempty(uTraj)
       JEffort = 0;
   else
       JEffort = trapz((0:numel(uTraj)-1)*dt, (uTraj).^2);
   end

   if isinf(chargingTime)
       JTime = 3000;
   else
       JTime = chargingTime;
   end
   effortW = 1e-3; 
   J = JTime + effortW * JEffort;

   VTPeak = -inf; if ~isempty(VTraj), VTPeak = max(VTraj); end
   TPeak  = -inf; if ~isempty(TTraj),  TPeak  = max(TTraj);  end
   IPeak  = 0;    if ~isempty(uTraj),  IPeak  = max(abs(uTraj)); end

   if isempty(VTraj)
       adherenceFrac = 1.0;
   else
       withinV = (VTraj <= mpcOpt.VT_max) & (VTraj >= mpcOpt.VT_min);

       TSample = interp1(data.t, TTraj, data.t(1:length(VTraj)), 'previous', 'extrap');
       withinT = (TSample <= mpcOpt.T_max);

       if isempty(uTraj)
           withinI = true(size(withinV));
       else
           u_sample = interp1(data.t(1:end-1), uTraj, data.t(1:length(VTraj)), 'previous', 'extrap');
           withinI = abs(u_sample) <= mpcOpt.I_max;
       end

       adherenceFrac = mean(withinV & withinT & withinI);
   end

   if numel(uTraj) < 3
       smoothnessU = 0;
   else
       smoothnessU = sum(abs(diff(uTraj)));
   end

   SOCEnd = SOCTraj(end);
   GSOC   = SOCEnd - SOCTarget;

   metrics = struct('chargingTime', chargingTime, ...
                    'adherenceFrac', adherenceFrac, ...
                    'JTime', JTime, ...
                    'JEffort', JEffort, ...
                    'smoothnessU', smoothnessU, ...
                    'VTPeak', VTPeak, ...
                    'TPeak', TPeak, ...
                    'IPeak', IPeak, ...
                    'SOCEnd', SOCEnd, ...
                    'GSOC', GSOC);

   if ~isempty(VTraj) && ~isempty(TTraj)
       phiSample = rbfFeatures([VTraj(1); TTraj(1)], mpcOpt.rbf);
       fprintf('diag: ||phi||=%.3e, phiMean=%.3e, phiMax=%.3e\n', ...
               norm(phiSample), mean(phiSample), max(phiSample));
       fprintf('diag: thetaNorm=%.3f, thetaMin=%.3f, thetaMax=%.3f\n', ...
               norm(theta), min(theta), max(theta));
       fprintf('diag: sample lRBF = %.3e\n', theta(:).' * phiSample(:));
   else
       fprintf('diag: empty VT/T trajectory, skipping RBF diagnostics\n');
   end

   disp(['SOCEnd = ', num2str(SOCEnd), ' | SOCTarget = ', num2str(SOCTarget)]);
end
