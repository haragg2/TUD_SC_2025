function [A_ineq, b_ineq, I_x, I_x0, I_ref] = build_constraints(params, dim)
    b_ineq = zeros(33, 1);
    A = zeros(29, 12);
    
    % Constraint 1
    A(1,:) = [1, 0, 0, params.xmax - 50, 0, 0, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(1) = params.xmax;
    
    % Constraint 2
    A(2,:) = [-1, 0, 0, -(50 + 1e-8), 0, 0, 0, 0, 0, 0, 0, 0]; % x, delta2
    b_ineq(2) = -(1e-8 + 50);
    
    % Constraint 3
    A(3,:) = [1, 0, 0, 0, params.xmax - 100, 0, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(3) = params.xmax;
    
    % Constraint 4
    A(4,:) = [-1, 0, 0, 0, -(100 + 1e-8), 0, 0, 0, 0, 0, 0, 0]; % x, delta3
    b_ineq(4) = -(1e-8 + 100);
    
    % Constraint 5
    A(5,:) = [1, 0, 0, 0, 0, params.xmax - 200, 0, 0, 0, 0, 0, 0]; % x
    b_ineq(5) = params.xmax;
    
    % Constraint 6
    A(6,:) = [-1, 0, 0, 0, 0, -(200 + 1e-8), 0, 0, 0, 0, 0, 0]; % x, delta3
    b_ineq(6) = -(1e-8 + 200);
    
    % Constraint 7
    A(7,:) = [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0]; % x, delta5
    b_ineq(7) = 0;

    % Constraint 8
    A(8,:) = [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0]; % x, delta5
    b_ineq(8) = 0;
    
    % Constraint 9
    A(9,:) = [0, 1, params.vmax - params.vg + 1e-8, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % v, delta1
    b_ineq(9) = params.vmax;
    
    % Constraint 10
    A(10,:) = [0, -1, -params.vg, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % v, delta1
    b_ineq(10) = -params.vg;
    
    % Constraint 11
    A(11,:) = [0, 1, 0, 0, 0, 0, params.vmax - params.alpha, 0, 0, 0, 0, 0]; % v, delta6
    b_ineq(11) = params.vmax;
    
    % Constraint 12
    A(12,:) = [0, -1, 0, 0, 0, 0, -params.alpha - 1e-8, 0, 0, 0, 0, 0]; % v, delta6
    b_ineq(12) = -(1e-8 + params.alpha);
    
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

    % Constraint 27 (delta_1 + delta_5 <= 1)
    A(27,:) = [0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0];
    b_ineq(27) = 0;
    
    % Constraint 28 (reference constraint - 1)
    % A(27,:) = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0];
    A(28,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0];
    b_ineq(28) = 0;
    
    % Constraint 29 (reference constraint - 2)
    % A(28,:) =  [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0];
    A(29,:) =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0];
    b_ineq(29) = 0;
    
    % % Constraint 30 (acceleration constraint - 1)
    b_ineq(30) = params.Ts * (params.acc_comf);

    % % Constraint 31 (acceleration constraint - 2)
    b_ineq(31) = params.Ts * (params.acc_comf);

    % % Constraint 32 (aux var for u constraint - 1)
    b_ineq(32) = 0;

    % % Constraint 33 (aux var for u constraint - 2)
    b_ineq(33) = 0;
    
   
    % Constraint matrices for k = 1 to Np
    I_delta = [A(1:29, 3:7); zeros(4, 5)];
    I_z = [A(1:29, 8:9); zeros(4, 2)];
    I_p = [A(1:29, 10); zeros(2, 1); -1; -1];
    I_q = [A(1:29, 11); zeros(4, 1)];

    I_delta = kron(eye(dim.Np), I_delta);
    I_z = kron(eye(dim.Np), I_z);
    I_p = kron(eye(dim.Np), I_p);
    I_q = kron(eye(dim.Np), I_q);

    I_x1 = [A(1:29, 1:2); [0, -1; 0, 1]; zeros(2, 2)];
    I_x2 = [zeros(27, 2); [0, 1; 0, -1; 0, 1; 0, -1]; zeros(2, 2)];

    I_x1p = kron(eye(dim.Np+1), I_x1);
    I_x1p = I_x1p(1:end-33, 3:end);

    I_x = kron(eye(dim.Np), I_x2) + I_x1p; % multiply this with x(1) ... x(dim.Np)

    I_x0 = [I_x1; zeros((dim.Np - 1) * size(I_x1, 1), 2)];

    I_u1 = [A(1:29, 12); zeros(2, 1); -1; 1];
    I_u2 = [zeros(31, 1); 1; -1];

    I_u2p = kron(eye(dim.Np+1), I_u2);
    I_u2p = I_u2p(1:end-33, 2:end);

    I_u = kron(eye(dim.Np), I_u1) + I_u2p % multiply this with u(0) ... u(dim.Np-1)

    I_ref = [zeros(27, 1); 1; -1; zeros(4, 1)]; % multiply this with vref
    I_ref = kron(eye(dim.Np), I_ref);

    b_ineq = repmat(b_ineq, dim.Np, 1);
    A_ineq = [I_z I_delta I_u I_q I_p];

% --------------------Inequality constraints ------------------------------%
% I_delta * [delta_1(1) ... delta_6(1) ... delta_1(dim.Np) ... delta_6(dim.Np)]^T +
% I_z * [z_1(1) ... z_2(1) ... z_1(dim.Np) ... z_2(dim.Np)]^T + 
% I_p * [p(1) ... p(dim.Np)]^T + I_q * [q(1) ... q(dim.Np)]^T + 
% I_u * [u(0) ... u(dim.Np-1)]^T + I_x * [x(1) ... x(dim.Np)]^T + 
% <= b_ineq + I_ref * [vref(1) ... vref(dim.Np)]^T - I_x0 * x(0);
%
% [I_z I_delata I_u I_q I_p] * [zp deltap up qp pp]^T <= -I_x * [x(1) ...
% x(dim.Np)]^T + b_ineq + I_ref * [vref(1) ... vref(dim.Np)]^T - I_x0 * x(0);
%
% [I_z I_delata I_u I_q I_p] * [zp deltap up qp pp]^T <= -I_x * (xp+)^T 
% + b_ineq + I_ref * [vref(1) ... vref(dim.Np)]^T - I_x0 * x(0);
%
% (xp+) = Ap * xp + Bp * [zp deltap up]^T + fp
%
% A_ineq * x <= b_ineq - I_x * (Ap * x(0) + Bp * x + fp) + I_ref * vref - I_x0 * x(0);
% (A_ineq + I_x * Bp) * x <= (b_ineq - I_x * fp) - (I_x * Ap + I_x0) * x(0) + I_ref * vref;
% ------------------------------------------------------------------------%

end
