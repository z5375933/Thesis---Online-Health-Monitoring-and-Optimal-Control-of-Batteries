paramsPlant = batteryParams();
params = batteryParams();

params.R0 = @(z) 1.7 * paramsPlant.R0(z);
params.R1 = @(z) 0.5 * paramsPlant.R1(z);
params.C1 = @(z) 1.3 * paramsPlant.C1(z);

x0List = {
    [0.2; 0; 298];
    [0.4; 0; 305]
    % [0.3; 0; 310]
};


mpcOpt = struct( ...
  'N', 5, ...
  'gamma_T', 1e3, ...
  'gamma_VT', 1e3, ...
  'T_max', 318, ...
  'VT_max', 4.2, ...
  'VT_min', 2.5, ...
  'I_max', 6.0, ...
  'I_min', 0.0, ...
  'dT', 2.0, ...
  'dVT', 0.05, ...
  'gamma_I', 1e-3, ...
  'I_nominal', 2.0, ...
  'gamma_dI', 1e-4, ...
  'rbf', buildCenters() ...
);

simOpt = struct('z_target', 0.8, 'M', 2000);

% BO options
BOOpt = struct();
BOOpt.N_iter = 15;
BOOpt.Ncand = 2000;
BOOpt.kappa = 1.5;
BOOpt.beta  = 1.3;
BOOpt.retrain_every = 1;

nTheta = size(mpcOpt.rbf.centers, 2);
BOOpt.theta_lb = -1 * ones(nTheta, 1);
BOOpt.theta_ub =  1 * ones(nTheta, 1);
BOOpt.fallback_scale = 0.03 * (BOOpt.theta_ub - BOOpt.theta_lb);

results = safeBOLoop(x0List, paramsPlant, params, mpcOpt, simOpt, BOOpt);

safeInds = find(results.safe_idx);
if isempty(safeInds)
    warning('No safe points found – falling back to global best J.');
    [~, best_idx] = min(results.JLog);
else
    [~, relBest] = min(results.JLog(safeInds));
    best_idx = safeInds(relBest);
end