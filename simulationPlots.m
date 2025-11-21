%% Sanity checks
assert(exist('paramsPlant','var')==1, 'paramsPlant not found – run Safe BO script first.');
assert(exist('params','var')==1,       'params not found – run Safe BO script first.');
assert(exist('mpcOpt','var')==1,      'mpcOpt not found – run Safe BO script first.');
assert(exist('simOpt','var')==1,      'simOpt not found – run Safe BO script first.');
assert(exist('x0List','var')==1,      'x0List not found – run Safe BO script first.');
assert(exist('thetaFinal','var')==1,  'thetaFinal not found – run Safe BO script first.');
assert(exist('results','var')==1,      'results not found – run Safe BO script first.');


% RL things (only needed for RL plots)
haveRL = exist('Q','var')==1 && exist('SOCBins','var')==1 && ...
         exist('TBins','var')==1 && exist('actionValues','var')==1 && ...
         exist('paramsEnv','var')==1 && exist('VTMax','var')==1 && ...
         exist('VTMin','var')==1 && exist('TMax','var')==1 && ...
         exist('IMax','var')==1 && exist('zTarget','var')==1 && ...
         exist('dt','var')==1;

%% 1) Baseline MPC closed-loop trajectories (IC1 & IC2)
disp('Plotting baseline MPC trajectories');

thetaZero = zeros(size(thetaFinal));
simOptBase = simOpt; simOptBase.zTarget = simOpt.zTarget;

for ic = 1:numel(x0List)
    x0 = x0List{ic};
    dataBase = simulateClosedLoop(x0, thetaZero, paramsPlant, params, mpcOpt, simOptBase);
    fig = figure;
    plotClosedLoop(dataBase, sprintf('Baseline MPC (IC %d)', ic));
    sgtitle(sprintf('Baseline MPC closed-loop (IC %d)', ic));
end


%% 2) BO-optimised MPC trajectories (IC1 & IC2)
disp('Plotting BO-optimised MPC trajectories');

for ic = 1:numel(x0List)
    x0 = x0List{ic};
    dataBO = simulateClosedLoop(x0, thetaFinal, paramsPlant, params, mpcOpt, simOpt);
    fig = figure;
    plotClosedLoop(dataBO, sprintf('BO–MPC (IC %d)', ic));
    sgtitle(sprintf('BO–MPC closed-loop (IC %d)', ic));
end


%% 3) RBF weight evolution across BO iterations
disp('Plotting RBF weight evolution');


thetaLg = results.thetaLog;      % [nIter x nTheta]
[nIter, nTheta] = size(thetaLg);


figure;
plot(1:nIter, thetaLg, 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('\theta_i');
title('RBF weight evolution');
grid on;


%% 4) Cost breakdown: Baseline MPC vs BO–MPC for IC1
disp('Plotting cost breakdown (IC1)');


x0 = x0List{1};


% Baseline
[~,~,~,~,~, runsBase] = thetaMultiIC(thetaZero, {x0}, ...
    paramsPlant, params, mpcOpt, simOpt);
mBase = runsBase{1}.metrics;


% BO–MPC
[~,~,~,~,~, runsBO] = thetaMultiIC(thetaFinal, {x0}, ...
    paramsPlant, params, mpcOpt, simOpt);
mBO = runsBO{1}.metrics;


JTime     = [mBase.JTime, mBO.JTime];
JEffort   = [mBase.JEffort, mBO.JEffort];
smoothness = [mBase.smoothnessU, mBO.smoothnessU];

figure;

subplot(3,1,1);
bar(JTime);
set(gca,'XTickLabel',{'Baseline','BO–MPC'});
ylabel('J_{time} (s)','Interpreter','tex');
title('Cost breakdown (IC1)');


subplot(3,1,2);
bar(JEffort);
set(gca,'XTickLabel',{'Baseline','BO–MPC'});
ylabel('J_{effort}','Interpreter','tex');


subplot(3,1,3);
bar(smoothness);
set(gca,'XTickLabel',{'Baseline','BO–MPC'});
ylabel('SmoothnessU','Interpreter','tex');
xlabel('Controller');


%% 5) Charging time comparison: Baseline vs BO–MPC (IC1 & IC2)
disp('Plotting charging time comparison (MPC vs BO–MPC)');

tBase = zeros(numel(x0List),1);
tBO   = zeros(numel(x0List),1);

for ic = 1:numel(x0List)
    x0 = x0List{ic};
    [~,~,~,~,~, runsBase] = thetaMultiIC(thetaZero, {x0}, paramsPlant, params, mpcOpt, simOpt);
    [~,~,~,~,~, runsBO]   = thetaMultiIC(thetaFinal, {x0}, paramsPlant, params, mpcOpt, simOpt);
    tBase(ic) = runsBase{1}.metrics.chargingTime;
    tBO(ic)   = runsBO{1}.metrics.chargingTime;
end

figure;
bar([tBase, tBO]);
set(gca,'XTickLabel',{'IC1','IC2'});
ylabel('Charging time (s)');
legend('Baseline MPC','BO–MPC','Location','northwest');
title('Charging time comparison (MPC vs BO–MPC)');


%% 6) RL plots (if RL info exists)
if haveRL
    disp('RL variables found – plotting RL-related figures');

    plotRLPolicyHeatmap(Q, SOCBins, TBins, actionValues);

    maxSteps = maxStepsPerEpisode; 
    rlTimes = zeros(numel(x0List),1);
    boTimes = zeros(numel(x0List),1);

    for ic = 1:numel(x0List)
        x0 = x0List{ic};

        rlMetrics = runRLEpisode(x0, Q, SOCBins, TBins, ...
                                    paramsEnv, VTMax, VTMin, ...
                                    TMax, IMax, zTarget, maxSteps, actionValues);

        [~,~,~,~,~, runsBO] = thetaMultiIC(thetaFinal, {x0}, paramsPlant, params, mpcOpt, simOpt);
        mBO = runsBO{1}.metrics;
        dataBO = runsBO{1}.data;


        rlTimes(ic) = rlMetrics.chargingTime;
        boTimes(ic) = mBO.chargingTime;


        plotRLvsMPCOverlay(rlMetrics, dataBO, sprintf('IC %d',ic));
    end

    figure;
    bar([boTimes, rlTimes]);
    set(gca,'XTickLabel',{'IC1','IC2'});
    ylabel('Charging time (s)');
    legend('BO–MPC','RL','Location','northwest');
    title('Charging time comparison (BO–MPC vs RL)');

else
    warning('RL variables not found');
end

disp('All plotting routines finished.');


%% Helper functions
function plotClosedLoop(data, subtitle_str)
    t  = data.t;
    z  = data.xTrj(1,:);
    T  = data.xTrj(3,:);
    I  = data.uTrj;
    VT = data.VT;


    subplot(4,1,1);
    plot(t, z, 'LineWidth', 1.5);
    ylabel('SOC');
    grid on;


    subplot(4,1,2);
    plot(t, T, 'LineWidth', 1.5);
    ylabel('T (K)');
    grid on;


    subplot(4,1,3);
    stairs(t(1:end-1), I, 'LineWidth', 1.5);
    ylabel('I (A)');
    grid on;


    subplot(4,1,4);
    plot(t(1:end-1), VT, 'LineWidth', 1.5);
    ylabel('V_T (V)');
    xlabel('Time (s)');
    grid on;


    if nargin > 1
        sgtitle(subtitle_str);
    end
end


function plotRLPolicyHeatmap(Q, SOCBins, TBins, actionValues)
    nZ = numel(SOCBins)-1;
    nT = numel(TBins)-1;
    [~, policy_idx] = max(Q,[],2); % best action per state
    policyMap = reshape(policy_idx, [nZ, nT])'; % T rows x Z cols


    figure;
    imagesc(SOCBins(1:end-1), TBins(1:end-1), policyMap);
    set(gca,'YDir','normal');
    colorbar;
    xlabel('SOC bin');
    ylabel('T bin (K)');
    title('RL learned policy (action index)');
end


function plotRLvsMPCOverlay(rlMetrics, dataBO, labelStr)
    tMpc  = dataBO.t;
    zMpc  = dataBO.x_trj(1,:);
    TMpc  = dataBO.x_trj(3,:);
    IMpc  = dataBO.u_trj;
    VTMpc = dataBO.VT;


    tRl  = rlMetrics.t;
    zRl  = rlMetrics.zHist;
    TRl  = rlMetrics.THist;
    IRl  = rlMetrics.IHist;
    VTRl = rlMetrics.VTHist;


    figure;
    subplot(4,1,1);
    plot(tMpc, zMpc, 'LineWidth', 1.5); hold on;
    plot(tRl, zRl, '--', 'LineWidth', 1.5);
    ylabel('SOC'); grid on;
    legend('BO–MPC','RL','Location','southeast');
    title(['Closed-loop comparison: ', labelStr]);


    subplot(4,1,2);
    plot(tMpc, TMpc, 'LineWidth', 1.5); hold on;
    plot(tRl, TRl, '--', 'LineWidth', 1.5);
    ylabel('T (K)'); grid on;


    subplot(4,1,3);
    stairs(tMpc(1:end-1), IMpc, 'LineWidth', 1.5); hold on;
    stairs(tRl(1:end-1), IRl, '--', 'LineWidth', 1.5);
    ylabel('I (A)'); grid on;


    subplot(4,1,4);
    plot(tMpc(1:end-1), VTMpc, 'LineWidth', 1.5); hold on;
    plot(tRl, VTRl, '--', 'LineWidth', 1.5);
    ylabel('V_T (V)'); xlabel('Time (s)'); grid on;
end
