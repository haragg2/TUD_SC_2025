clc;
clear;
close all;
addpath('functions/');

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% Model Parameters
m_w = 1.5;          % kg; mass of wheel
m_p = 6.16;         % kg; mass of pendulum
l = 0.46;           % m; distance between COG of pendulum and center of wheel
r = 0.054;          % m; radius of wheel
g = 9.81;           % m/s^2; gravity
J = 0.1;            % kgm^2; wheel centroidal inertia?

% The slope can be modified
beta = 10 * (pi/180);       % 10 degrees; slope angle

% Desired Reference Parameters (for non-liniear system)
theta_eq = 10;              % rad; wheel displacement
theta_d_eq = 0;
theta_p_eq = asin((m_w + m_p) * r * sin(beta) / (m_p * l));    % 1.73 degrees; pendulum desired angle
theta_pd_eq = 0;

Ts = 0.05;           % s; sampling time
T_sim = 10;          % s; simulation time
T = 0:Ts:T_sim;      % simulation time steps

% Linearizing the model
sys_cont = sys_dyn(m_w, m_p, l, r, g, beta, J);

% Discretize the system
sys_dis = c2d(sys_cont, Ts, 'zoh');

LTI.A = sys_dis.A;
LTI.B = sys_dis.B;
LTI.C = sys_dis.C;
LTI.D = sys_dis.D;

%% Bounds

xlb = [-30; (-1 /(r)); (-8*pi/180); (-pi/10)];
xub = [30; (1 /(r)); (8*pi/180); (pi/10)];

ulb = [-2];
uub = [2];

N = 10; % Horizon
%% LQR

% Tuning weights
Q = 1 * eye(size(LTI.A));            % state cost
R = 16 * eye(length(LTI.B(1,:)));    % input cost

% Controllability check
ctrb_sys = ctrb(LTI.A, LTI.B);
unco = length(LTI.A) - rank(ctrb_sys);
if unco > 0
    warning('Discretized linear model is uncontrollable');
end

% Observability check
obs_sys = obsv(LTI.A, LTI.C);
unob = length(LTI.A) - rank(obs_sys);
if unob > 0
    warning('Discretized linear model is unobservable');
end

% Check if (A, Q') has no unobservable modes on the unit circle
isNoUnobservableModes = rank(ctrb(LTI.A, transpose(Q))) == size(LTI.A,1);
if isNoUnobservableModes == 0
    warning("Pair (A, Q') has some unobservable modes");
end


[S_lqr_cost, K_lqr, ~] = idare(LTI.A, LTI.B, Q, R);
% [K_lqr, S_lqr_cost] = dlqr(LTI.A, LTI.B, Q, R);
K_lqr = -K_lqr;
A_K = LTI.A + LTI.B * K_lqr;

% x0 = [-27.2; 17.7; -0.13; 0.21]; % In terminal set Xf
x0 = [29.4; -18.1; 0.13; -0.28]; % In terminal set Xf

sysd_lqr = ss(A_K, [], LTI.C, LTI.D, Ts);
u = zeros(size(T));
lsim(sysd_lqr, u, T, x0);

[~, t_LQ, x_LQ] = lsim(sysd_lqr, u, T, x0);
u_LQ = K_lqr * x_LQ';

figure;
plot(t_LQ, u_LQ);
title('Control Input for LQ Control');
xlabel('Time (seconds)');
ylabel('Control Input');
grid on;

%% Compute X_f.

Xn = struct();
V = struct();
Z = struct();

[Xn.('lqr'), V.('lqr'), Z.('lqr')] = findXn(LTI.A, LTI.B, K_lqr, N, xlb, xub, ulb, uub, 'lqr');

Xf = Xn.lqr{1};
X1 = Xn.lqr{2};
X2 = Xn.lqr{10};
%% Plot Xf and Xn

% Generate 100000 random points within the specified bounds
num_points = 1000000;
x_samples = bsxfun(@plus, xlb, bsxfun(@times, rand(numel(xlb), num_points), (xub - xlb)));

% Check if Ax <= b is satisfied for each point
satisfied_points_Xf = all(Xf.A * x_samples <= Xf.b, 1);
satisfied_points_X1 = all(X1.A * x_samples <= X1.b, 1);
satisfied_points_X2 = all(X2.A * x_samples <= X2.b, 1);

% Plot the points that satisfy the condition on the x1-x2 plane
figure;
hold on;
plot(x_samples(1, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(1, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(1, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$\theta_{w}$');
ylabel('$\theta_{p}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;


figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{w}}$$');
ylabel('$$\theta_{p}$$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(4, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{p}}$$');
ylabel('$\theta_{p}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(1, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(1, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(1, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{w}}$$');
ylabel('$\theta_{w}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;


%% Linear MPC - Regulation

% Solves the following optimization problem:
%
%      N-1
% min  sum [ 1/2(u_k' H u_k + d u_k ) ]
% u_k  k=0
%
%    s.t    G x_k + H u_k <= psi (states and control input constraints)
%           A_f x_N <= b_f (terminal constraint)

Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));

% Defines the dimensions
dim.N = N;                      % Horizon
dim.nx = size(LTI.A, 1);        % Number of states
dim.nu = size(LTI.B, 2);        % Number of inputs

% Verified - its correct
[P, S, W] = predmodgen(LTI,dim);

x0 = [16.9, 3.89, -0.11, -0.09]';                       % initial state for MPC (should be from set Xn)
% x0 = [29.4; -18.1; 0.13; -0.28]; % In terminal set Xf

H = (S'*Qbar*S + Rbar + 2*W'*S_lqr_cost*W);             % Hessian for quadratic cost in inputs
d = (x0'*P'*Qbar*S + 2*x0'*(LTI.A^N)'*S_lqr_cost*W)';   % linear term for cost in inputs

G_N = kron(eye(N), Z.('lqr').('G'));
H_N = kron(eye(N), Z.('lqr').('H'));
psi_N = repmat(Z.('lqr').('psi'), 10, 1);

x = zeros(length(LTI.A(:,1)), size(T, 2));      % state trajectory
x(:,1) = x0;                                    % set first state as x0

u = zeros(length(LTI.B(1,:)), size(T, 2));      % control inputs
% y = zeros(length(LTI.C(:,1)),T);              % measurements 
t = zeros(1, size(T, 2));                       % time vector
LV = zeros(1, size(T, 2));                      % Lyapunov function value at every x
Vf = zeros(1, size(T, 2));
% l = zeros(1, size(T, 2));

for k = 1:1:size(T, 2)
    t(k) = (k-1) * Ts;
    if ( mod(t(k), 1) == 0 )
        fprintf('t = %d sec \n', t(k));
    end

    % determine reference states based on reference input r
    x0 = x(:, k);
    d = (x0'*P'*Qbar*S + 2*x0'*(LTI.A^N)'*S_lqr_cost*W)';

    % % compute control action
    % cvx_begin quiet
    %     variable u_N(1*N)
    %     minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
    % 
    %     % state-input constraints
    %     (G_N * S + H_N) * u_N <= psi_N - G_N * (P*x0);
    %     % terminal constraints
    %     Xf.A * W * u_N <= Xf.b - Xf.A * LTI.A^N * x0;
    % cvx_end

    Options = sdpsettings('verbose', 0, 'solver', 'quadprog');
    u_N_ = sdpvar(N, 1);                % define optimization variable

    Constraint = [(G_N * S + H_N) * u_N_ <= (psi_N - G_N * (P*x0));
                  Xf.A * W * u_N_ <= Xf.b - Xf.A * LTI.A^N * x0];   %define constraints

    Objective = 0.5*(u_N_' * H * u_N_) + d' * u_N_;  %define cost function

    diagnostic = optimize(Constraint, Objective, Options);  %solve the problem
    u_N = value(u_N_);                  %assign the solution to uopt
    FVAL = 0.5 * (u_N' * H * u_N) + d' * u_N;
    status.optimal = diagnostic.problem;

    u(:,k) = u_N(1:1); % MPC control action

    % apply control action
    x(:,k+1) = LTI.A*x(:,k) + LTI.B*u(:,k);
    % y(:,k) = C*x(:,k);

    LV(k) = FVAL;
    Vf(k) = 0.5 * (x(:,k)' * S_lqr_cost * x(:,k));
    % l(k) = 0.5*x(:,k)'*Q*x(:,k) + 0.5*u(:,k)'*R*u(:,k);
end

%%  Simulation of regulation
close all;
x0 = [16.9, 3.89, -0.11, -0.09]';
figure;
lsim(sysd_lqr, u, T, x0);

xw = @(tt) (r * cos(beta) * [1 0 0 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1)));
yw = @(tt) (r * sin(beta) * [1 0 0 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1)));

xp = @(tt) ((r * cos(beta) * [1 0 0 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1))) + (l * sin([0 0 1 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1)))));
yp = @(tt) ((r * sin(beta) * [1 0 0 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1))) + (l * cos([0 0 1 0]*(LTI.A* x(:,int8(tt/0.05)+1) + LTI.B*u(int8(tt/0.05)+1)))));

figure;
axis equal;
hold on;
fanimator(@(tt) plot(xp(tt), yp(tt),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(1) xw(tt)],[yw(1) yw(tt)],'b-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(tt) xp(tt)],[yw(tt) yp(tt)],'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot(xw(tt), yw(tt),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) text(-0.3,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',20);
hold off;

%% Function definitions

function sys_cont = sys_dyn(m_w, m_p, l, R, g, beta, J)
 
    theta_pd = asin((m_w + m_p) * R * sin(beta) / (m_p * l));    % 1.73 degrees; pendulum desired angle at eq
     
    a = J + ((m_w + m_p) * R^2);
    b_eq = m_p * R * l * cos(theta_pd + beta);
    c =  m_p * l^2;
    %d = (m_w + m_p) * g * R * sin(beta);    % At eq, tau = d
     
    A_1 = m_p * g * l / ((a * c) - b_eq^2);
     
    A = [0 1 0 0;
        0 0 (-b_eq) * A_1 0; %0 0 (-b_eq - c) * A_1 0;
        0 0 0 1;
        0 0 (a) * A_1 0]; %0 0 (a - b_eq) * A_1 0];
     
    B = [0;
        (c + b_eq) / ((a * c) - b_eq^2);
        0;
        (-a - b_eq) / ((a * c) - b_eq^2)];
     
    C = eye(4);
     
    D = 0;
     
    sys_cont = ss(A,B,C,[]);
end


function [P,S,W] = predmodgen(LTI,dim)

    % Prediction matrices generation
    % This function computes the prediction matrices to be used in the
    % optimization problem

    % Prediction matrix from initial state
    P = zeros(dim.nx * dim.N, dim.nx);
    for k = 0:dim.N-1
        P(k*dim.nx+1:(k+1)*dim.nx,:) = LTI.A^k;
    end

    % Prediction matrix from input
    S = zeros(dim.nx * dim.N, dim.nu * dim.N);
    for k = 1:dim.N-1
        for i = 0:k-1
            S(k*dim.nx+1:(k+1)*dim.nx, i*dim.nu+1:(i+1)*dim.nu) = LTI.A^(k-1-i) * LTI.B;
        end
    end

    W = zeros(dim.nx, dim.nu * dim.N);
    for i = 0:dim.N-1
        W(:, i*dim.nu+1:(i+1)*dim.nu) = LTI.A^(dim.N-i-1) * LTI.B;
    end

end