function results = safeBOLoop(x0List, paramsPlant, params, mpcOpt, simOpt, BOOpt)

rng('shuffle');

thetaLog = [];    % Log of parameter vectors
JLog     = [];    % Log of cost values 
GiLog    = [];    % Log of safety margins [G1, G2, G3]
safeIdx  = [];    % Log of safe evaluations
dataLog  = {};    
maxIter = BOOpt.N_iter;

fprintf('\n SAFE BO START \n');

theta0 = zeros(size(mpcOpt.rbf.centers, 2), 1);
nSeed = 8;
seed_scale = 0.03;
nTheta = numel(theta0);

[J0, VT0, T0, I0, G4_0, runs0] = thetaMultiIC(theta0, x0List, paramsPlant, params, mpcOpt, simOpt);
G1_0 = mpcOpt.VT_max - VT0;
G2_0 = mpcOpt.T_max  - T0;
G3_0 = mpcOpt.I_max  - I0;

thetaLog = [thetaLog; theta0.'];
JLog     = [JLog; J0];
GiLog    = [GiLog; [G1_0, G2_0, G3_0, G4_0]];
safeIdx  = [safeIdx; (G1_0 >= 0) && (G2_0 >= 0) && (G3_0 >= 0) && (G4_0 >= 0)];
dataLog{end+1} = runs0;

for si = 1:nSeed
   th = theta0 + seed_scale * randn(nTheta,1);

   [J_s, VT_s, T_s, I_s, G4_s, runs_s] = ...
       thetaMultiIC(th, x0List, paramsPlant, params, mpcOpt, simOpt);

   G1_s = mpcOpt.VT_max - VT_s;
   G2_s = mpcOpt.T_max  - T_s;
   G3_s = mpcOpt.I_max  - I_s;

   thetaLog = [thetaLog; th.'];
   JLog     = [JLog; J_s];
   GiLog    = [GiLog; [G1_s, G2_s, G3_s, G4_s]];
   safeIdx  = [safeIdx; (G1_s >= 0) && (G2_s >= 0) && (G3_s >= 0) && (G4_s >= 0)];
   dataLog{end+1} = runs_s; 
end

gpJ  = fitrgp(thetaLog, JLog, 'KernelFunction','ardsquaredexponential', ...
             'InitialSigma',0.02, ...
             'SigmaLowerBound', 1e-4, ...
             'Standardize',true);

gpG1 = fitrgp(thetaLog, GiLog(:,1), 'KernelFunction','ardsquaredexponential', ...
             'InitialSigma',0.02, ...
             'SigmaLowerBound', 1e-4, ...
             'Standardize',true);

gpG2 = fitrgp(thetaLog, GiLog(:,2), 'KernelFunction','ardsquaredexponential', ...
             'InitialSigma',0.02, ...
             'SigmaLowerBound', 1e-4, ...
             'Standardize',true);

gpG3 = fitrgp(thetaLog, GiLog(:,3), 'KernelFunction','ardsquaredexponential', ...
             'InitialSigma',0.02, ...
             'SigmaLowerBound', 1e-4, ...
             'Standardize',true);

gpG4 = fitrgp(thetaLog, GiLog(:,4), 'KernelFunction','ardsquaredexponential', ...
             'InitialSigma',0.02, ...
             'SigmaLowerBound', 1e-4, ...
             'Standardize',true);

% BO options with defaults
Ncand = BOOpt.Ncand;
kappa = BOOpt.kappa;
beta  = BOOpt.beta;
retrain_every = BOOpt.retrain_every;
lb = BOOpt.theta_lb(:).';
ub = BOOpt.theta_ub(:).';
fallback_scale = BOOpt.fallback_scale(:).';

for iter = 1:maxIter
   fprintf('\nIteration %d / %d\n', iter, maxIter);

   if iter <= 3
       thetaNew = 0.1 * randn(size(theta0));
       thetaNew = min(max(thetaNew, lb.'), ub.');
       fprintf('Random exploration (cold start)\n');
   else
       % Safe BO: propose new theta using GP predictions
       thetaNew = proposeNextTheta(gpJ, {gpG1, gpG2, gpG3, gpG4}, ...
           thetaLog, JLog, GiLog, ...
           struct('Ncand',Ncand, 'kappa',kappa, 'beta',beta, ...
                  'lb',lb, 'ub',ub, 'fallback_scale',fallback_scale));
       thetaNew = min(max(thetaNew(:), lb.'), ub.');
   end

    [J_new, VT_new_peak, T_new_peak, I_new_peak, G4_new, runsNew] = ...
        thetaMultiIC(thetaNew, x0List, paramsPlant, params, mpcOpt, simOpt);
    
    G1New = mpcOpt.VT_max - VT_new_peak;
    G2New = mpcOpt.T_max  - T_new_peak;
    G3New = mpcOpt.I_max  - I_new_peak;
    
    isSafe = (G1New >= 0) && (G2New >= 0) && (G3New >= 0) && (G4_new >= 0);
    
    thetaLog = [thetaLog; thetaNew.'];
    JLog     = [JLog; J_new];
    GiLog    = [GiLog; G1New, G2New, G3New, G4_new];
    safeIdx  = [safeIdx; isSafe];
    dataLog{end+1} = runsNew;  

   if mod(iter, retrain_every) == 0
       gpJ  = fitrgp(thetaLog, JLog, 'KernelFunction','ardsquaredexponential', ...
                     'InitialSigma',0.02, ...
                     'SigmaLowerBound', 1e-4, ...
                     'Standardize',true);

       gpG1 = fitrgp(thetaLog, GiLog(:,1), 'KernelFunction','ardsquaredexponential', ...
                     'InitialSigma',0.02, ...
                     'SigmaLowerBound', 1e-4, ...
                     'Standardize',true);

       gpG2 = fitrgp(thetaLog, GiLog(:,2), 'KernelFunction','ardsquaredexponential', ...
                     'InitialSigma',0.02, ...
                     'SigmaLowerBound', 1e-4, ...
                     'Standardize',true);

       gpG3 = fitrgp(thetaLog, GiLog(:,3), 'KernelFunction','ardsquaredexponential', ...
                     'InitialSigma',0.02, ...
                     'SigmaLowerBound', 1e-4, ...
                     'Standardize',true);

       gpG4 = fitrgp(thetaLog, GiLog(:,4), 'KernelFunction','ardsquaredexponential', ...
                     'InitialSigma',0.02, ...
                     'SigmaLowerBound', 1e-4, ...
                     'Standardize',true);
   end

    fprintf('‣ |θ|=%.3f | J=%.6f | VT=%.3f | T=%.3f | I=%.3f | Safe=%d\n', ...
            norm(thetaNew), J_new, VT_new_peak, T_new_peak, I_new_peak, isSafe);
    if size(thetaLog,1) > 1
        prevTheta = thetaLog(end-1,:).';
        fprintf('Δθ norm = %.6f\n', norm(thetaNew - prevTheta));
    end
end

[~, best_idx] = min(JLog);
best_runs    = dataLog{best_idx};  
best_metrics = best_runs{1}.metrics; 

fprintf('\n BEST SAFE BO RESULT \n');
fprintf('Iteration %d | Best J=%.3f\n', best_idx, JLog(best_idx));
fprintf('Charging time (IC1) = %.2f s\n', best_metrics.charging_time);
fprintf('Adherence (IC1) = %.3f\n', best_metrics.adherence_frac);
fprintf('VT_peak (IC1) = %.3f | T_peak (IC1) = %.3f | I_peak (IC1) = %.3f\n', ...
    best_metrics.VT_peak, best_metrics.T_peak, best_metrics.I_peak);
fprintf('J_time (IC1) = %.3f | J_effort (IC1) = %.3f | Smoothness_u (IC1) = %.3f\n', ...
    best_metrics.J_time, best_metrics.J_effort, best_metrics.smoothness_u);

results = struct('theta_log', thetaLog, 'J_log', JLog, 'Gi_log', GiLog, ...
                 'safe_idx', safeIdx, 'gpJ', gpJ, 'gpG1', gpG1, 'gpG2', gpG2, 'gpG3', gpG3, 'gpG4', gpG4, ...
                 'data', {dataLog}, 'options', BOOpt);


chargingTimes = zeros(length(results.data),1);
for i = 1:length(results.data)
   runs_i = results.data{i};           % cell array over ICs
   metrics_i1 = runs_i{1}.metrics;     % IC1
   if isfield(metrics_i1,'charging_time') && ~isnan(metrics_i1.charging_time)
       chargingTimes(i) = metrics_i1.charging_time;
   else
       chargingTimes(i) = inf;
   end
   fprintf('Charging time at iteration %d (IC1): %.2f s\n', i, chargingTimes(i));
end

fprintf('\n SAFE BO COMPLETE \n');


