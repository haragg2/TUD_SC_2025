% clear all
close all
% clc
addpath('functions/');
addpath('linear_system_matrices/');

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Linearize the system

% global MPC_data;

% Linearization point
theta1_eq = pi;
theta2_eq = 0;

x_eq = [theta1_eq; 0; theta2_eq; 0];

sys_mat = load("LinearSys_0_0.mat");
sys_ct = ss(sys_mat.sys_A, sys_mat.sys_B, [1, 0, 0, 0; 0, 0, 1, 0], []);

%% Discretize the system

h = 0.01;           % s; sampling time
T_sim = 3;          % s; simulation time
T = 0:h:T_sim;      % simulation time steps

sys_dis = c2d(sys_ct, h);

% System matrices
A = sys_dis.A;
B = sys_dis.B;
C = sys_dis.C;
D = sys_dis.D;

%% Bounds

xlb = [-0.3; -100; -0.4; -100];
xub = -xlb;
ulb = -1;
uub = 1;

N = 5; % Horizon

% Defines the dimensions
dim.nx = size(sys_dis.A, 1);        % Number of states
dim.nu = size(sys_dis.B, 2);        % Number of inputs
dim.ny = size(sys_dis.C, 1);        % Number of outputs

%% LQR

% Tuning weights
Q = [40 0 0 0;
     0 10 0 0;
     0 0 70 0;
     0 0 0 10];

R = 0.1;

% Find LQR gain matrix
[K, P] = dlqr(A, B, Q, R);
K = -K; % Sign convention

% Controllability check
ctrb_sys = ctrb(A, B);
unco = length(A) - rank(ctrb_sys);
if unco > 0
    warning('Discretized linear model is uncontrollable');
end

% Observability check
obs_sys = obsv(A, C);
unob = length(A) - rank(obs_sys);
if unob > 0
    warning('Discretized linear model is unobservable');
end

% Check if (A, Q') has no unobservable modes on the unit circle
isNoUnobservableModes = rank(ctrb(A, transpose(Q))) == size(A,1);
if isNoUnobservableModes == 0
    warning("Pair (A, Q') has some unobservable modes");
end

A_K = A + B * K;

x0 = [0.15, 0, -0.1, 0]';       % initial condition for Xn
u = zeros(size(T));                     % 0 input
sysd_lqr = ss(A_K, [], C, D, h);       

figure;
grid on;
lsim(sysd_lqr, u, T, x0);

[~, t_LQ, x_LQ] = lsim(sysd_lqr, u, T, x0);
u_LQ = K * x_LQ';

figure;
stairs(t_LQ, u_LQ);
title('Control Input for LQ Control');
xlabel('Time (seconds)');
ylabel('Control Input');
grid on;

%% Compute X_f

Xn = struct();
V = struct();
Z = struct();

[Xn.('lqr'), V.('lqr'), Z.('lqr')] = findXn(A, B, K, N, xlb, xub, ulb, uub, 'lqr');

%% Plot Xf and Xn
Xf = Xn.lqr{1};
X1 = Xn.lqr{floor(N/2)};
X2 = Xn.lqr{N+1};

% % Generate 500000 random points within the specified bounds
% num_points = 500000;
% x_samples = bsxfun(@plus, xlb, bsxfun(@times, rand(numel(xlb), num_points), (xub - xlb)));
% % x_samples = [x_samples, [0.2654; -3.0552; -0.1543; 5.2605]];
% 
% % Check if Ax <= b is satisfied for each point
% satisfied_points_Xf = all(Xf.A * x_samples <= Xf.b, 1);
% satisfied_points_X1 = all(X1.A * x_samples <= X1.b, 1);
% satisfied_points_X2 = all(X2.A * x_samples <= X2.b, 1);

% % Plot the points that satisfy the condition on the 2D plane
% figure;
% sgtitle('Projection of points satisfying $X.A*x <= X.b$ on 2D plane');
% hplots = gobjects(3, 1);
% subplot(2, 2, 1);
% hold on;
% hplots(1) = plot(x_samples(1, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'ro', 'MarkerSize', 3);
% hplots(2) = plot(x_samples(1, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
% hplots(3) = plot(x_samples(1, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
% hold off;
% xlabel('$\theta_{1}$');
% ylabel('$\theta_{2}$');
% grid on;
% 
% subplot(2, 2, 2);
% hold on;
% hplots(1) = plot(x_samples(2, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'ro', 'MarkerSize', 3);
% hplots(2) = plot(x_samples(2, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
% hplots(3) = plot(x_samples(2, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
% hold off;
% xlabel('$$\dot{\theta_{1}}$$');
% ylabel('$$\theta_{2}$$');
% grid on;
% 
% subplot(2, 2, 3);
% hold on;
% hplots(1) = plot(x_samples(4, satisfied_points_X2), x_samples(1, satisfied_points_X2), 'ro', 'MarkerSize', 3);
% hplots(2) = plot(x_samples(4, satisfied_points_X1), x_samples(1, satisfied_points_X1), 'go', 'MarkerSize', 3);
% hplots(3) = plot(x_samples(4, satisfied_points_Xf), x_samples(1, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
% hold off;
% xlabel('$$\dot{\theta_{2}}$$');
% ylabel('$\theta_{1}$');
% grid on;
% 
% subplot(2, 2, 4);
% hold on;
% hplots(1) = plot(x_samples(2, satisfied_points_X2), x_samples(1, satisfied_points_X2), 'ro', 'MarkerSize', 3);
% hplots(2) = plot(x_samples(2, satisfied_points_X1), x_samples(1, satisfied_points_X1), 'go', 'MarkerSize', 3);
% hplots(3) = plot(x_samples(2, satisfied_points_Xf), x_samples(1, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
% hold off;
% xlabel('$$\dot{\theta_{1}}$$');
% ylabel('$\theta_{1}$');
% grid on;
% hL = legend(hplots, {'$X_{15}$', '$X_{7}$', '$X_f$'}, 'Interpreter', 'latex', 'Orientation', 'horizontal');
% 
% % Set the location of the legend to 'southoutside' which positions
% % it below the subplots and centers it.
% set(hL, 'Location', 'southoutside', 'Box', 'off');

%% Regulation MPC

model_mpc = struct('A', A, 'B', B, 'C', C, 'Bd', zeros(size(B)), 'Cd', zeros(size(C, 1), 1), 'N', N);
constraint = Z.lqr;
penalty = struct('Q', Q, 'R', R, 'P', P);
terminal = Xn.lqr{1}; % LQR terminal set

% x0 = [0.25, 0, -0.1, 0]';  % initial condition for XN
x0 = [0.2654; -3.0552; -0.1543; 5.2605];
xr = [0; 0; 0; 0]; % reference x_r set to 0 for regulation
ref = [repmat([xr;0],N,1); xr];

x = zeros(dim.nx, size(T, 2));
x(:,1) = x0;

usim = zeros(dim.nu, size(T, 2)-1);
V_N = zeros(size(T, 2)-1,1);
V_f = zeros(size(T, 2)-1,1);
l_xu = zeros(size(T, 2)-1,1);
l_xu0 = zeros(size(T, 2)-1,1);
t = zeros(size(T, 2)-1,1);
mpcmats = []; % Calculated first time and then reused

%% for actual implementation on simulink

MPC_data = struct('model', model_mpc, 'constraint', constraint, 'penalty', penalty, 'terminal', terminal, 'mpcmats', mpcmats);
save('MPC_DATA.mat', "MPC_data");
%% Functions

% function [xr, ur] = targetSelector(LTI, Z, dim, d_hat, yref)
% 
%     eqconstraints.A = [eye(dim.nx) - LTI.A, -LTI.B; LTI.C, zeros(size(LTI.C, 1), dim.nu)];
%     eqconstraints.b = [LTI.Bd * d_hat; yref - (LTI.Cd*d_hat)];
% 
%     ineqconstraints.A = [Z.('G'), Z.('H')];
%     ineqconstraints.b = Z.('psi');
% 
%     H = blkdiag(zeros(dim.nx), eye(dim.nu));
%     h = zeros(dim.nx+dim.nu, 1);
% 
%     options1 = optimoptions(@quadprog);
%     options1.OptimalityTolerance=1e-20;
%     options1.ConstraintTolerance=1.0000e-15;
%     options1.Display='off';
%     xur=quadprog(H,h,ineqconstraints.A,ineqconstraints.b,eqconstraints.A,eqconstraints.b,[],[],[],options1);
%     xr = xur(1:dim.nx);
%     ur = xur(dim.nx+1:end);
% end
