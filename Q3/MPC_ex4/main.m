clear all
close all
clc

%% Model definition

% Model Parameters

m_w = 1.5;          % kg; mass of wheel
m_p = 6.16;         % kg; mass of pendulum
l = 0.46;           % m; distance between COG of pendulum and center of wheel
r = 0.054;          % m; radius of wheel
g = 9.81;           % m/s^2; gravity
J = 0.1;            % kgm^2; wheel centroidal inertia

% The slope can be modified
beta = 10 * (pi/180);       % 10 degrees; slope angle

% Desired Reference Parameters (for non-linear system)
theta_eq = 10;              % rad; wheel displacement
theta_d_eq = 0;
theta_p_eq = asin((m_w + m_p) * r * sin(beta) / (m_p * l));    % pendulum desired angle
theta_pd_eq = 0;

Ts = 0.05;           % s; sampling time
T_sim = 10;          % s; simulation time
T = 0:Ts:T_sim;      % simulation time steps

% Linearizing the model
sys_cont = sys_dyn(m_w, m_p, l, r, g, beta, J);

% Discretize the system
sys_dis = c2d(sys_cont, Ts, 'zoh');

A = sys_dis.A;
B = sys_dis.B;
C = sys_dis.C;
D = sys_dis.D;

%% Bounds

xlb = [-30; (-1 /(r)); (-8*pi/180); (-pi/10)];
xub = [30; (1 /(r)); (8*pi/180); (pi/10)];

ulb = [-2];
uub = [2];

N = 5; % Horizon

% Defines the dimensions
dim.nx = size(sys_dis.A, 1);        % Number of states
dim.nu = size(sys_dis.B, 2);        % Number of inputs
dim.ny = size(sys_dis.C, 1);        % Number of outputs

%% LQR

% Tuning weights
Q = 1 * eye(size(A));            % state cost
R = 16 * eye(length(B(1,:)));    % input cost

% Find LQR

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

x0 = [29.4; -18.1; 0.13; -0.28]; % In terminal set Xf
%x0 = [16.9, 3.89, -0.11, -0.09]'; % in Set XN - control input exceeds constraints
sysd_lqr = ss(A_K, [], C, D, Ts);
u = zeros(size(T));
figure;
grid on;
lsim(sysd_lqr, u, T, x0);

[~, t_LQ, x_LQ] = lsim(sysd_lqr, u, T, x0);
u_LQ = K * x_LQ';

figure;
plot(t_LQ, u_LQ);
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

% Generate 100000 random points within the specified bounds
num_points = 100000;
x_samples = bsxfun(@plus, xlb, bsxfun(@times, rand(numel(xlb), num_points), (xub - xlb)));

% Check if Ax <= b is satisfied for each point
satisfied_points_Xf = all(Xf.A * x_samples <= Xf.b, 1);
satisfied_points_X1 = all(X1.A * x_samples <= X1.b, 1);
satisfied_points_X2 = all(X2.A * x_samples <= X2.b, 1);

% Plot the points that satisfy the condition on the x1-x2 plane
figure;
hold on;
plot(x_samples(1, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'ro', 'MarkerSize', 3);
plot(x_samples(1, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(1, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$\theta_{w}$');
ylabel('$\theta_{p}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'ro', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{w}}$$');
ylabel('$$\theta_{p}$$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(4, satisfied_points_X2), x_samples(3, satisfied_points_X2), 'ro', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_X1), x_samples(3, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(4, satisfied_points_Xf), x_samples(3, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{p}}$$');
ylabel('$\theta_{p}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;

figure;
hold on;
plot(x_samples(2, satisfied_points_X2), x_samples(1, satisfied_points_X2), 'ro', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_X1), x_samples(1, satisfied_points_X1), 'go', 'MarkerSize', 3);
plot(x_samples(2, satisfied_points_Xf), x_samples(1, satisfied_points_Xf), 'bo', 'MarkerSize', 3);
hold off;
xlabel('$$\dot{\theta_{w}}$$');
ylabel('$\theta_{w}$');
title('Points satisfying Ax <= b on x1-x2 plane');
grid on;


%% Regulation MPC

model = struct('A', A, 'B', B, 'N', N);
constraint = Z.lqr;
penalty = struct('Q', Q, 'R', R, 'P', P);
terminal = Xn.lqr{1}; % LQR terminal set

x0 = [16.9, 3.89, -0.11, -0.09]'; % in Set XN
xr = [0; 0; 0; 0]; % referece x set to 0

x = zeros(dim.nx, size(T, 2));
x(:,1) = x0;

usim = zeros(dim.nu, size(T, 2)-1);
V_N = zeros(size(T, 2)-1,1);
mpcmats = []; % Calculated first time and then reused


for t = 1:1:size(T, 2)-1
    model.x0 = x(:,t);
    [xk, uk, FVAL, status, mpcmats] = linearmpc(xr,0,0,model, constraint, penalty, ...
                                             terminal, mpcmats);
    usim(:,t) = uk(:,1);
    V_N(t) = FVAL;
    tspan = [T(t), T(t+1)];
    x(:,t+1) = NLSystemDynamics(xk(:,1), tspan, usim(:,t));
end

figure;
plot(V_N);
title("Optimal Cost funcion V_N^0");
grid on;

figure;
title("Actual State Response");
subplot(2, 2, 1);
plot(x(1,:), 'LineWidth', 2), grid on;
subplot(2, 2, 2);
plot(x(2,:), 'LineWidth', 2), grid on;
subplot(2, 2, 3);
plot(x(3,:), 'LineWidth', 2), grid on;
subplot(2, 2, 4);
plot(x(4,:), 'LineWidth', 2), grid on;

figure;
plot(usim, 'LineWidth', 2), grid on;
title("Control Input");

%%  Simulation of regulation

% Non-linear equations
xw = @(tt) (r * cos(beta) * x(1, (uint16(tt/0.05)+1)));
yw = @(tt) (r * sin(beta) * x(1, (uint16(tt/0.05)+1)));

xp = @(tt) ((r * cos(beta) * x(1, (uint16(tt/0.05)+1)) + (l * sin(x(3, (uint16(tt/0.05)+1))))));
yp = @(tt) ((r * sin(beta) * x(1, (uint16(tt/0.05)+1)) + (l * cos(x(3, (uint16(tt/0.05)+1))))));

scaling_factor = 0.6;
angle = -pi/3; % Example angle of 45 degrees

figure;
axis equal;
hold on;
fanimator(@(tt) plot(xp(tt), yp(tt),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(Ts) xw(tt)],[yw(Ts) yw(tt)],'b-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot(xw(tt), yw(tt),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% Plot the L-shaped line
fanimator(@(tt) plot([xw(tt) xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle)], [yw(tt) yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle) xp(tt)], [yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle) yp(tt)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% fanimator(@(tt) text(-0.3,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',20);
hold off;


%% MPC - Reference Tracking

yref = [6; 0; 0; 0];

% Calculate xref, uref
eqconstraints.A = [eye(dim.nx) - A, -B; C, zeros(dim.ny, dim.nu)];
eqconstraints.b = [zeros(dim.nx, 1); yref];

ineqconstraints.A = [Z.('lqr').('G'), Z.('lqr').('H')];
ineqconstraints.b = Z.('lqr').('psi');

H = blkdiag(zeros(dim.nx), eye(dim.nu));
h = zeros(dim.nx+dim.nu, 1);

options1 = optimoptions(@quadprog); 
options1.OptimalityTolerance=1e-20;
options1.ConstraintTolerance=1.0000e-15;
options1.Display='iter';
xur=quadprog(H,h,ineqconstraints.A,ineqconstraints.b,eqconstraints.A,eqconstraints.b,[],[],[],options1);
xr = xur(1:dim.nx);
ur = xur(dim.nx+1:end);
ref = [repmat([xr;ur],N,1); xr];

x0 = [16.9, 3.89, -0.11, -0.09]'; % in Set XN
x(:,1) = x0; % initial condition
usim = zeros(dim.nu, size(T, 2)-1);
mpcmats = []; % Calculated first time and then reused

for t = 1:1:size(T, 2)-1
    model.x0 = x(:,t);
    [xk, uk, FVAL, status, mpcmats] = linearmpc(xr,ref,0,model, constraint, penalty, ...
                                             terminal, mpcmats);
    usim(:,t) = uk(:,1);
    V_N(t) = FVAL;
    tspan = [T(t), T(t+1)];
    x(:,t+1) = NLSystemDynamics(xk(:,1), tspan, usim(:,t));
end

figure;
plot(V_N);
title("Optimal Cost funcion V_N^0");
grid on;

figure;
title("Actual State Response");
subplot(2, 2, 1);
plot(x(1,:), 'LineWidth', 2), grid on;
subplot(2, 2, 2);
plot(x(2,:), 'LineWidth', 2), grid on;
subplot(2, 2, 3);
plot(x(3,:), 'LineWidth', 2), grid on;
subplot(2, 2, 4);
plot(x(4,:), 'LineWidth', 2), grid on;

figure;
plot(usim, 'LineWidth', 2), grid on;
title("Control Input");
grid on;

%% Simulation of reference

% Non-linear equations
xw = @(tt) (r * cos(beta) * x(1, (uint16(tt/0.05)+1)));
yw = @(tt) (r * sin(beta) * x(1, (uint16(tt/0.05)+1)));

xp = @(tt) ((r * cos(beta) * x(1, (uint16(tt/0.05)+1)) + (l * sin(x(3, (uint16(tt/0.05)+1))))));
yp = @(tt) ((r * sin(beta) * x(1, (uint16(tt/0.05)+1)) + (l * cos(x(3, (uint16(tt/0.05)+1))))));

figure;
axis equal;
hold on;
fanimator(@(tt) plot(xp(tt), yp(tt),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(Ts) xw(tt)],[yw(Ts) yw(tt)],'b-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot(xw(tt), yw(tt),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% Plot the L-shaped line
fanimator(@(tt) plot([xw(tt) xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle)], [yw(tt) yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle) xp(tt)], [yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle) yp(tt)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% fanimator(@(tt) text(-0.3,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',20);
hold off;

%% Output feedback

Bd = B;
C_out = [1 0 0 0; 0 0 1 0; 0 0 0 1];
Cd = [0.5; 0; 0];

if ((rank([eye(dim.nx) - A -Bd; C_out Cd]) - dim.nx - 1) ~= 0) || (rank(obsv(A,C_out)) ~= dim.nx)
    warning ('Augmented state is not observable');
end

aug_sys.A = [A Bd; zeros(1,dim.nx) eye(1,size(Bd,2))];
aug_sys.B = [B;zeros(1,size(B,2))];
aug_sys.C = [C_out Cd];

obs_poles = [0.5*pole(sysd_lqr); 0.5];
L_obs = place(aug_sys.A', aug_sys.C', obs_poles)';

if abs(eig(aug_sys.A - L_obs*aug_sys.C)) >= 1
    warning ('Observer is not stable');
end

disturbance = 1.5;
d_hat = 0;

yref = [10; 0; 0];

xr = zeros(dim.nx,1);
ur = [];

ye = zeros(length(yref),size(T, 2)-1);
ye(:,1)=aug_sys.C*[x0; disturbance];

xehat = zeros(dim.nx+1, size(T, 2));
xehat(:,1)=[0, 0, 0, 0, d_hat];

x0 = [1, 0, 0.08, 0]';
x(:,1) = x0;

usim = zeros(dim.nu, size(T, 2)-1);
d_est = zeros(1,size(T, 2)-1);
mpcmats = [];

H = blkdiag(zeros(dim.nx), eye(dim.nu));
h = zeros(dim.nx+dim.nu, 1);
    
options1 = optimoptions(@quadprog); 
options1.OptimalityTolerance=1e-20;
options1.ConstraintTolerance=1.0000e-15;
options1.Display='off';

for t = 1:1:size(T, 2)-1

    d_hat = xehat(end,t);
    d_est(t) = d_hat;
    dist = repmat((Bd*d_hat),N,1);

    eqconstraints.A = [eye(dim.nx) - A, -B; C_out, zeros(size(C_out, 1), dim.nu)];
    eqconstraints.b = [Bd*d_hat; yref-(Cd*d_hat)];
    
    xur=quadprog(H,h,[],[],eqconstraints.A,eqconstraints.b,[],[],[],options1);
    xr = xur(1:dim.nx);
    ur = xur(dim.nx+1:end);
    ref = [repmat([xr;ur],N,1); xr];
    
    model.x0 = x(:,t); % xehat(1:end-1,t);

    [xk, uk, FVAL, status] = linearmpc(xr,ref, dist, model, constraint, penalty, ...
                                             terminal, []);
    usim(:,t) = uk(:,1);

    tspan = [T(t), T(t+1)];
    x(:,t+1) = NLSystemDynamics(xk(:,1), tspan, usim(:,t)+disturbance); % given that Bd = B
    ye(:,t) = aug_sys.C * [x(:,t); disturbance];
    xehat(:,t+1) = aug_sys.A*xehat(:,t) + aug_sys.B*usim(:,t) + L_obs*(ye(:,t) - (aug_sys.C*xehat(:,t)));
    V_N(t) = FVAL;
   
end

figure;
plot(d_est);
title("Estimated Disturbance");
grid on;

figure;
plot(V_N);
title("Optimal Cost funcion V_N^0");
grid on;

figure;
title("Estimated State Response");
subplot(2, 2, 1);
plot(xehat(1,:), 'LineWidth', 2), grid on;
subplot(2, 2, 2);
plot(xehat(2,:), 'LineWidth', 2), grid on;
subplot(2, 2, 3);
plot(xehat(3,:), 'LineWidth', 2), grid on;
subplot(2, 2, 4);
plot(xehat(4,:), 'LineWidth', 2), grid on;

figure;
title("Actual State Response");
subplot(2, 2, 1);
plot(x(1,:), 'LineWidth', 2), grid on;
subplot(2, 2, 2);
plot(x(2,:), 'LineWidth', 2), grid on;
subplot(2, 2, 3);
plot(x(3,:), 'LineWidth', 2), grid on;
subplot(2, 2, 4);
plot(x(4,:), 'LineWidth', 2), grid on;


figure;
plot(ye(1,:), 'LineWidth', 2), grid on;
title("Output");

figure;
plot(usim, 'LineWidth', 2), grid on;
title("Control Input");

%% Simulation of output feedback

% Non-linear equations
xw = @(tt) (r * cos(beta) * x(1, (uint16(tt/0.05)+1)));
yw = @(tt) (r * sin(beta) * x(1, (uint16(tt/0.05)+1)));

xp = @(tt) ((r * cos(beta) * x(1, (uint16(tt/0.05)+1)) + (l * sin(x(3, (uint16(tt/0.05)+1))))));
yp = @(tt) ((r * sin(beta) * x(1, (uint16(tt/0.05)+1)) + (l * cos(x(3, (uint16(tt/0.05)+1))))));

figure;
axis equal;
hold on;
fanimator(@(tt) plot(xp(tt), yp(tt),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(Ts) xw(tt)],[yw(Ts) yw(tt)],'b-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot(xw(tt), yw(tt),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% Plot the L-shaped line
fanimator(@(tt) plot([xw(tt) xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle)], [yw(tt) yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(tt)-scaling_factor*l*cos(x(3, (uint16(tt/0.05)+1))-angle) xp(tt)], [yw(tt)+scaling_factor*l*sin(x(3, (uint16(tt/0.05)+1))-angle) yp(tt)] ,'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);

% fanimator(@(tt) text(-0.3,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',20);
hold off;

%% Functions

function sys_cont = sys_dyn(m_w, m_p, l, R, g, beta, J)
 
    theta_pd = asin((m_w + m_p) * R * sin(beta) / (m_p * l));    % pendulum desired angle at eq
     
    a = J + ((m_w + m_p) * R^2);
    b_eq = m_p * R * l * cos(theta_pd + beta);
    c =  m_p * l^2;
    %d = (m_w + m_p) * g * R * sin(beta);    % At eq, tau = d
     
    A_1 = m_p * g * l / ((a * c) - b_eq^2);
     
    A = [0 1 0 0;
        0 0 -b_eq * A_1 0; %0 0 (-b_eq - c) * A_1 0;
        0 0 0 1;
        0 0 a * A_1 0]; %0 0 (a - b_eq) * A_1 0];
     
    B = [0;
        (c + b_eq) / ((a * c) - b_eq^2);
        0;
        (-a - b_eq) / ((a * c) - b_eq^2)];
     
    C = eye(4);
     
    D = 0;
     
    sys_cont = ss(A,B,C,D);
end

function [X] = NLSystemDynamics(x0, tspan, u)
    syms theta(t) theta_p(t) m_w m_p g l R B J Tau

    theta_d = diff(theta);
    theta_pd = diff(theta_p);

    theta_dd = diff(theta_d);
    theta_pdd = diff(theta_pd);

    eq1 = theta_dd(t) == (((m_p * l^2) + (m_p * R * l * cos(theta_p(t) + B))) * Tau ...
        + ((m_p * l^2) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2) ...
        - (m_p * R * l * cos(theta_p(t) + B) * m_p * g * l * sin(theta_p(t))) ...
        - (m_p * l^2) * ((m_w + m_p) * g * R * sin(B))) / (((J ...
        + ((m_w + m_p) * R^2)) * (m_p * l^2)) - ((m_p * R * l * cos(theta_p(t) + B))^2));
    
    
    eq2 = theta_pdd(t) == (((-(J + ((m_w + m_p) * R^2)) - (m_p * R * l * cos(theta_p(t) + B))) * Tau) ...
        + ((J + ((m_w + m_p) * R^2)) * m_p * g * l * sin(theta_p(t))) ...
        - ((m_p * R * l * cos(theta_p(t) + B)) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2) ...
        + ((m_p * R * l * cos(theta_p(t) + B)) * ((m_w + m_p) * g * R * sin(B))))  / (((J ...
        + ((m_w + m_p) * R^2)) * (m_p * l^2)) - ((m_p * R * l * cos(theta_p(t) + B))^2));
    
    % Model Parameters
    m_w = 1.5;          % kg; mass of wheel
    m_p = 6.16;         % kg; mass of pendulum
    l = 0.46;           % m; distance between COG of pendulum and center of wheel
    R = 0.054;          % m; radius of wheel
    g = 9.81;           % m/s^2; gravity
    J = 0.1;            % kgm^2; wheel centroidal inertia
    B = 10 * (pi/180);

    eq1 = subs(eq1);
    eq2 = subs(eq2);

    [V, ~] = odeToVectorField(eq2, eq1);
    M = matlabFunction(V,'vars',{'t','Y', 'Tau'});
    
    x0 = x0';
    theta_p_eq = asin((m_w + m_p) * R * sin(B) / (m_p * l));
    x0(3) = x0(3) + theta_p_eq;  % Add the eqbm pendulum angle to NL system 
    
    Tau = u  + (m_w + m_p) * g * R * sin(B);
    [~, x_step] = ode45(@(t, Y)M(t, Y, Tau), tspan, x0);
    X = x_step(end, :)';
    X(3) = X(3) - theta_p_eq; % Subtract the eqbm pendulum angle to feed it back to the linear system
end