clc;
clear;
addpath('functions/');

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

% Linearizing the model
sys_cont = sys_dyn(m_w, m_p, l, r, g, beta, J);

% Discretize the system
sys_dis = c2d(sys_cont, Ts, 'zoh');

% Controllability check
ctrb_sys = ctrb(sys_dis.A, sys_dis.B);
unco = length(sys_dis.A) - rank(ctrb_sys);
if unco > 0
    warning('Discretized linear model is uncontrollable');
end

% Observability check
obs_sys = obsv(sys_dis.A, sys_dis.C);
unob = length(sys_dis.A) - rank(obs_sys);
if unob > 0
    warning('Discretized linear model is unobservable');
end

LTI.A = sys_dis.A;
LTI.B = sys_dis.B;
LTI.C = sys_dis.C;

% Tuning weights
Q = 1 * eye(size(LTI.A));            % state cost
R = 1 * eye(length(LTI.B(1,:)));    % input cost

% [S_lqr_cost, K_lqr, ~] = idare(LTI.A, LTI.B, Q, R);
[K_lqr, S_lqr_cost] = dlqr(LTI.A, LTI.B, Q, R);
K_lqr = -K_lqr;
A_K = LTI.A + LTI.B * K_lqr;


%% LQR



%%

% Bounds.
xlb = [-30; (-1 /(r)); (-8*pi/180); (-pi/10)];
xub = [30; (1 /(r)); (8*pi/180); (pi/10)];
% xlb = [-inf(); -inf(); -inf(); -inf()];
% xub = [inf(); inf(); inf(); inf()];
ulb = [-1];
uub = [1];
% ulb = [-inf()];
% uub = [inf()];

N = 10; % Horizon
%% Compute X_f.

Xn = struct();
V = struct();
Z = struct();

[Xn.('lqr'), V.('lqr'), Z.('lqr')] = findXn(LTI.A, LTI.B, K_lqr, N, xlb, xub, ulb, uub, 'lqr');
%%
Xf = Xn.lqr{1};
X1 = Xn.lqr{2};
X2 = Xn.lqr{9};

%%

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
xlabel('x1');
ylabel('x2');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;


figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('x1');
ylabel('x2');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(4, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('x1');
ylabel('x2');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(1, satisfied_points_X2), 'yo', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(1, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(1, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('x1');
ylabel('x2');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;





%%
% % computes a control invariant set for LTI system x^+ = A*x+B*u
% system = LTISystem('A', A_K);
% system.x.min = [-50; -50; -50; -50];
% system.x.max = [50; 50; 50; 50];
% system.u.min = [-10];
% system.u.max = [10];
% InvSet = system.invariantSet()
% InvSet.plot()



Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));

% Defines the dimensions
dim.N = N;
dim.nx = size(LTI.A, 1);
dim.nu = size(LTI.B, 2);
[P, Z, W] = predmodgen(LTI,dim);

x0 = [0.1, 0, 0.1, 0]';
H = (Z'*Qbar*Z + Rbar + 2*W'*S_lqr_cost*W);
d = (x0'*P'*Qbar*Z + 2*x0'*(LTI.A^N)'*S_lqr_cost*W)';

T = N/Ts;
x_lim_vec_full = repmat([5; 5; 5; 5], [N 1]);

x = zeros(length(LTI.A(:,1)),T);    % state trajectory
x(:,1) = x0;
u = zeros(length(LTI.B(1,:)),T);    % control inputs
% y = zeros(length(LTI.C(:,1)),T);    % measurements 
t = zeros(1,T);                 % time vector

u_limit = 2;
% for k = 1:1:T
%     t(k) = (k-1) * Ts;
%     if ( mod(t(k), 1) == 0 )
%         fprintf('t = %d sec \n', t(k));
%     end
% 
%     % determine reference states based on reference input r
%     x0 = x(:,k);
%     d = (x0'*P'*Qbar*Z + 2*x0'*(LTI.A^N)'*S_lqr_cost*W)';
% 
%     % compute control action
%     cvx_begin quiet
%         variable u_N(1*N)
%         minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
%         % input constraints
%         u_N <=  u_limit*ones(1*N,1);
%         u_N >= -u_limit*ones(1*N,1);
%         % state constraints
%         Z*u_N <= -P*x0 + x_lim_vec_full;
%         Z*u_N >= -P*x0 - x_lim_vec_full; 
%     cvx_end
% 
%     u(:,k) = u_N(1:1); % MPC control action
% 
%     % apply control action
%     x(:,k+1) = LTI.A*x(:,k) + LTI.B*u(:,k); % + B_ref*r(:,k);
%     % y(:,k) = C*x(:,k);
% 
%     [X,K, ~] = idare(LTI.A,LTI.B,Q,R);
%     Vf(k) = 0.5*x(:,k)'*X*x(:,k);
%     l(k) = 0.5*x(:,k)'*Q*x(:,k) + 0.1*u(:,k)'*R*u(:,k);
% end
% 
% plot(t, x(1, 1:200));

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


function [P,Z,W] = predmodgen(LTI,dim)

    % Prediction matrices generation
    % This function computes the prediction matrices to be used in the
    % optimization problem

    % Prediction matrix from initial state
    P = zeros(dim.nx * dim.N, dim.nx);
    for k = 0:dim.N-1
        P(k*dim.nx+1:(k+1)*dim.nx,:) = LTI.A^k;
    end

    % Prediction matrix from input
    Z = zeros(dim.nx * dim.N, dim.nu * dim.N);
    for k = 1:dim.N-1
        for i = 0:k-1
            Z(k*dim.nx+1:(k+1)*dim.nx, i*dim.nu+1:(i+1)*dim.nu) = LTI.A^(k-1-i) * LTI.B;
        end
    end

    W = zeros(dim.nx, dim.nu * dim.N);
    for i = 0:dim.N-1
        W(:, i*dim.nu+1:(i+1)*dim.nu) = LTI.A^(dim.N-i-1) * LTI.B;
    end

end