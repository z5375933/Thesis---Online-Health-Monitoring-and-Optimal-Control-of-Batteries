function [JMean, VTPeakMax, TPeakMax, IPeakMax, G4Min, runs] = ...
   thetaMultiIC(theta, x0List, paramsPlant, params, mpcOpt, simOpt)

   nIC = numel(x0List);

   JVals       = zeros(nIC,1);
   VTPeaks     = zeros(nIC,1);
   TPeaks      = zeros(nIC,1);
   IPeaks      = zeros(nIC,1);
   G4Vals      = zeros(nIC,1);

   runs         = cell(nIC,1);

   for l = 1:nIC
       x0 = x0List{l};

       [J, VT_peak, T_peak, I_peak, data, metrics] = ...
           evaluateTheta(theta, x0, paramsPlant, params, mpcOpt, simOpt);

       JVals(l)   = J;
       VTPeaks(l) = VT_peak;
       TPeaks(l)  = T_peak;
       IPeaks(l)  = I_peak;
       G4Vals(l)  = metrics.G_SOC;

       runs{l} = struct('data', data, 'metrics', metrics);
   end

   JMean = mean(JVals);

   VTPeakMax = max(VTPeaks);
   TPeakMax  = max(TPeaks);
   IPeakMax  = max(IPeaks);

   G4Min = min(G4Vals);
end
