function xp = MLD_test(params, x0, u)
    
    
    %-------------------------System Dynamics (MLD)---------------------------%
    A1 = [1, params.Ts; 0, (1 - (params.Ts * params.Ps) / params.m)];
    
    B1 = [0, 0; 
        params.Ts * (params.r1 - params.r2), -(params.Ts / ...
        params.m) * ((params.beta/params.alpha) - params.Ps)];
    
    B2 = [zeros(1, 5);
        0, -(params.Ts * params.g * params.theta2), -( ...
        params.Ts * params.g * params.theta2), ( ...
        params.Ts * params.g * params.theta4), ...
        -(params.Ts / params.m) * (params.alpha * params.Ps - params.beta)];
    
    B3 = [0; params.Ts * params.r2];
    
    f = [0; (-params.Ts*params.g*params.theta4)-(params.Ts / params.m) * (params.beta - params.alpha * params.Ps)];
    
    
    %-------------------------------------------------------------------------%
    
    
    b_ineq = zeros(27, 1);
    A = zeros(27, 12);
    
    % Constraint 1
    A(1,:) = [1, 0, 0, params.xmax - 50, 0, 0, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(1) = params.xmax;
    
    % Constraint 2
    A(2,:) = [-1, 0, 0, -(50 + eps(1)), 0, 0, 0, 0, 0, 0, 0, 0]; % x, delta2
    b_ineq(2) = -(eps(1) + 50);
    
    % Constraint 3
    A(3,:) = [1, 0, 0, 0, params.xmax - 100, 0, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(3) = params.xmax;
    
    % Constraint 4
    A(4,:) = [-1, 0, 0, 0, -(100 + eps(1)), 0, 0, 0, 0, 0, 0, 0]; % x, delta3
    b_ineq(4) = -(eps(1) + 100);
    
    % Constraint 5
    A(5,:) = [1, 0, 0, 0, 0, params.xmax - 200, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(5) = params.xmax;
    
    % Constraint 6
    A(6,:) = [-1, 0, 0, 0, 0, -(200 + eps(1)), 0, 0, 0, 0, 0, 0]; % x, delta3
    b_ineq(6) = -(eps(1) + 200);
    
    % Constraint 7
    A(7,:) = [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0]; % x, delta5
    b_ineq(7) = 0;

    % Constraint 8
    A(8,:) = [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0]; % x, delta5
    b_ineq(8) = 0;
    
    % Constraint 9
    A(9,:) = [0, 1, params.vmax - params.vg + eps(1), 0, 0, 0, 0, 0, 0, 0, 0, 0]; % v, delta1
    b_ineq(9) = params.vmax;
    
    % Constraint 10
    A(10,:) = [0, -1, -params.vg, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % v, delta1
    b_ineq(10) = -params.vg;
    
    % Constraint 11
    A(11,:) = [0, 1, 0, 0, 0, 0, params.vmax - params.alpha, 0, 0, 0, 0, 0]; % v, delta6
    b_ineq(11) = params.vmax;
    
    % Constraint 12
    A(12,:) = [0, -1, 0, 0, 0, 0, -params.alpha - eps(1) + params.vg, 0, 0, 0, 0, 0]; % v, delta6
    b_ineq(12) = -(eps(1) + params.alpha);
    
    % Constraint 13
    A(13,:) = [0, 0, -params.umax, 0, 0, 0, 0, 1, 0, 0, 0, 0]; % z1, delta1
    b_ineq(13) = 0;
    
    % Constraint 14
    A(14,:) = [0, 0, params.umin, 0, 0, 0, 0, -1, 0, 0, 0, 0]; % z1, delta1
    b_ineq(14) = 0;
    
    % Constraint 15
    A(15,:) = [0, 0, -params.umin, 0, 0, 0, 0, 1, 0, 0, 0, -1]; % z1, delta1
    b_ineq(15) = -params.umin;
    
    % Constraint 16
    A(16,:) = [0, 0, params.umax, 0, 0, 0, 0, -1, 0, 0, 0, 1]; % z1, delta1
    b_ineq(16) = params.umax;
    
    % Constraint 17
    A(17,:) = [0, 0, 0, 0, 0, 0, -params.vmax, 0, 1, 0, 0, 0]; % z2
    b_ineq(17) = 0;
    
    % Constraint 18
    A(18,:) = [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0]; % z2
    b_ineq(18) = 0;
    
    % Constraint 19
    A(19,:) = [0, -1, 0, 0, 0,0, 0, 0, 1, 0, 0, 0]; % v
    b_ineq(19) = 0;
    
    % Constraint 20
    A(20,:) = [0, 1, 0, 0, 0, 0, params.vmax, 0, -1, 0, 0, 0];
    b_ineq(20) = params.vmax;
    
    % Constraint 21
    A(21,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
    b_ineq(21) = params.umax;
    
    % Constraint 22
    A(22,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];
    b_ineq(22) = -params.umin;
    
    % Constraint 23
    A(23,:) = [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    b_ineq(23) = 0;
    
    % Constraint 24
    A(24,:) = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    b_ineq(24) = params.vmax;
    
    % Constraint 25
    A(25,:) = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    b_ineq(25) = 0;
    
    % Constraint 26
    A(26,:) = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    b_ineq(26) = params.xmax;

    % Constraint 27
    A(27,:) = [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0];
    b_ineq(27) = 1;
    
    % Constraint 29 (acceleration constraint - 1)
    b_ineq(28) = params.Ts * (params.acc_comf);

    % Constraint 30 (acceleration constraint - 2)
    b_ineq(29) = params.Ts * (params.acc_comf);
    
    
    % Constraint matrices for k = 1 to Np
    
    % I_delta = [A(1:27, 3:7); zeros(2, 5)];
    % I_z = [A(1:27, 8:9); zeros(2, 2)];
    
    I_delta = A(1:27, 3:7);
    I_z = A(1:27, 8:9);
    
    % I_x1 = [A(1:27, 1:2); [0, -1; 0, 1]];
    I_x1 = A(1:27, 1:2);

    % I_x2 = [zeros(27, 2); [0, 1; 0, -1]];
    
    % I_u1 = [A(1:27, 12); zeros(2, 1)];
    I_u1 = A(1:27, 12);
    
    % Inequality constraint
    A_i = [I_z I_delta];
    b_i = b_ineq(1:27) - I_x1 * x0 - I_u1 * u;

    % A_i = [I_z I_delta] + I_x2 * [B1 B2];
    % b_i = b_ineq(1:29) - I_x2 * f - (I_x1 + I_x2*A1) * x0 - (I_u1 + I_x2 * B3) * u;
    
    
    options = optimoptions("intlinprog", "Display", "off");
    
    intcon = [3, 4, 5, 6, 7];
    lb = [-inf; -inf; zeros(5, 1)];
    ub = [inf; inf; ones(5, 1)];
    
    [output, fval, flag, msg] = intlinprog(zeros(1, 7), intcon, A_i, b_i, [], [], lb, ub, [], options);
    % disp(msg)
    
    if flag == 1
        z = output(1:2);
        delta = output(3:7);
        % x0
        xp = A1 * x0 + B1 * z + B2 * delta + B3 * u + f;
    else
        disp("Infeasible Solution. Previous value used");
        xp = x0;
    end

