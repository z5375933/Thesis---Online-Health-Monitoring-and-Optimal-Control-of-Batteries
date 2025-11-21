% clear; 
% clc; 
% close all;


%% 1. ENVIRONMENT SETUP

paramsEnv = batteryParams();
dt = paramsEnv.dt;

VTMax = 4.2;
VTMin = 2.5;
TMax  = 318;
IMax  = 6.0;
zTarget = 0.8;


x0List = {
    [0.2; 0; 298]
    [0.4; 0; 305]
    % [0.3; 0; 310]
};

maxStepsPerEpisode = 2000;

actionValues = [0;2;4;6];
nA = numel(actionValues);


%% 2. STATE DISCRETISATION

SOCBins = 0:0.05:1.0;
TBins   = 290:2:320;   

nZ = numel(SOCBins) - 1;
nT = numel(TBins)   - 1;
nS = nZ * nT;


%% 3. Q-LEARNING HYPERPARAMETERS


numEpisodes = 1500;
alpha = 0.1;
gamma = 0.99;
eps   = 1.0;
epsMin   = 0.05;
epsDecay = 0.995;

Q = zeros(nS, nA);
episodeRewards = zeros(numEpisodes,1);


%% 4. TRAINING LOOP

fprintf('\n Q-learning training start \n');


for ep = 1:numEpisodes
    x0 = rlReset(x0List);
    [~, VT0] = batteryStep(x0, 0, paramsEnv);
    s = [x0(1); VT0; x0(2); x0(3)];


    s_idx = discretizeState(s, SOCBins, TBins);
    done = false;
    totalReward = 0;


    for t = 1:maxStepsPerEpisode
        if rand < eps
            a_idx = randi(nA);
        else
            [~, a_idx] = max(Q(s_idx,:));
        end


        [s_next, r, done] = rlStep(s, a_idx, paramsEnv, ...
                                   VTMax, VTMin, TMax, IMax, ...
                                   zTarget, actionValues);


        totalReward = totalReward + r;


        sNext_idx = discretizeState(s_next, SOCBins, TBins);


        % Q-learning update
        Q(s_idx, a_idx) = (1 - alpha)*Q(s_idx, a_idx) + ...
                          alpha*(r + gamma*max(Q(sNext_idx,:)));


        s     = s_next;
        s_idx = sNext_idx;


        if done
            break;
        end
    end


    episodeRewards(ep) = totalReward;
    eps = max(epsMin, eps * epsDecay);


    if mod(ep,100)==0
        fprintf('Episode %d / %d, totalReward = %.1f, eps = %.3f, steps = %d\n', ...
            ep, numEpisodes, totalReward, eps, t);
    end
end


fprintf('\n Q-learning training complete \n');


%% 5. PLOT TRAINING REWARD


figure;
plot(episodeRewards,'LineWidth',1.5);
xlabel('Episode');
ylabel('Total reward');
title('Q-learning training progress');
grid on;


%% 6. EVALUATE RL POLICY VS MPC / BO-MPC (IF AVAILABLE)


fprintf('\n Evaluation: RL vs MPC \n');

haveMPC = exist('paramsPlant','var') && exist('params','var') && exist('mpcOpt','var') && exist('simOpt','var');


if haveMPC
    fprintf('MPC parameters found\n');
else
    fprintf('MPC parameters NOT found\n');
end


for ic = 1:numel(x0List)
    x0 = x0List{ic};

    rlMetrics = runRLEpisode(x0, Q, SOCBins, TBins, ...
                                paramsEnv, VTMax, VTMin, ...
                                TMax, IMax, zTarget, ...
                                maxStepsPerEpisode, actionValues);


    fprintf('\nIC %d (RL): time = %.1f s, violated = %d\n', ...
            ic, rlMetrics.chargingTime, rlMetrics.violated);

    figure;
    subplot(4,1,1);
    plot(rlMetrics.t, rlMetrics.zHist,'LineWidth',1.5);
    ylabel('SOC');
    title(sprintf('RL closed-loop trajectory (IC %d)', ic));


    subplot(4,1,2);
    plot(rlMetrics.t, rlMetrics.T_Hist,'LineWidth',1.5);
    ylabel('T (K)');


    subplot(4,1,3);
    stairs(rlMetrics.t(1:end-1), rlMetrics.IHist,'LineWidth',1.5);
    ylabel('I (A)');


    subplot(4,1,4);
    plot(rlMetrics.t, rlMetrics.VTHist,'LineWidth',1.5);
    ylabel('V_T (V)');
    xlabel('Time (s)');

    if haveMPC
        thetaBase = zeros(size(mpc_opt.rbf.centers,2),1);
        [~, ~, ~, ~, dataBase, metricsBase] = ...
            evaluateTheta(thetaBase, x0, paramsPlant, params, mpcOpt, simOpt);


        fprintf('IC %d (Baseline MPC): time = %.1f s, adherence = %.3f\n', ...
                ic, metricsBase.chargingTime, metricsBase.adherenceFrac);


        if exist('thetaFinal','var')
            [~, ~, ~, ~, dataBO, metricsBO] = ...
                evaluateTheta(thetaFinal, x0, paramsPlant, params, mpcOpt, simOpt);


            fprintf('IC %d (BO–MPC): time = %.1f s, adherence = %.3f\n', ...
                    ic, metricsBO.chargingTime, metricsBO.adherenceFrac);
        end
    end
end
