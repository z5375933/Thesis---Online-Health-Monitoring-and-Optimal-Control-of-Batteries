function [u_seq, x_pred, info] = mpcController(x0, theta, params, opt)

    N   = opt.N;
    nx  = 3;
    nu  = 1;
    dt  = params.dt;

    % Decision vector: [u_0 ... u_{N-1}]
    u0 = zeros(N,1);
    lb = opt.I_min * ones(N,1);
    ub = opt.I_max * ones(N,1);

    % fmincon options
    fopts = optimoptions('fmincon', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', false, ...
        'SpecifyConstraintGradient', false, ...
        'MaxFunctionEvaluations', 5e4, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-8);
    
    function J = obj(u)
        x = x0;
        J = 0;
        for k = 1:N
            Ik = u(k);
            [x, VT] = ecmDynamics(x, Ik, params);

            z = x(1);
            T = x(3);

            % Base MPC cost
            l_base = (0.8 - z)^2 ...
                   + opt.gamma_T  * max(0, T  - (opt.T_max  - opt.dT))^2 ...
                   + opt.gamma_VT * max(0, VT - (opt.VT_max - opt.dVT))^2;

            l_rbf = 0;
            if ~isempty(theta)
                phi = rbfFeatures([VT; T], opt.rbf);
                l_rbf = (theta(:).' * phi(:));
            end

            J = J + l_base + l_rbf;
        end
    end

    [u_star, fval, exitflag, output] = fmincon(@obj, u0, [], [], [], [], lb, ub, [], fopts);

    u_seq = u_star;
    info.exitflag = exitflag;
    info.output = output;
    info.cost = fval;

    % Predict trajectory
    x_pred = zeros(nx, N+1);
    x_pred(:,1) = x0;
    x = x0;
    for k = 1:N
        [x, ~] = ecmDynamics(x, u_seq(k), params);
        x_pred(:,k+1) = x;
    end
end