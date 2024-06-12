clear all;
close all;
clc;

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%%
syms m c v vmax alpha beta

p1 = (beta/alpha)*v;
p2 = ((c*vmax^2-beta) / (vmax-alpha))*(v - alpha) + beta;

a1_int = @(v) (c * v^2 - p1)^2;
a1 =  simplify(int(a1_int, v, 0, alpha));

a2_int = @(v) (c * v^2 - p2)^2;
a2 = simplify(int(a2_int, v, alpha, vmax));

pa = a1+a2;
dp_alpha = diff(pa, alpha);
dp_beta = diff(pa, beta);

% Solve the system of equations
solutions = solve([dp_alpha == 0, dp_beta == 0], [alpha, beta]);

alpha_sol = solutions.alpha;
beta_sol = solutions.beta;

% Parameters
g = 9.8;
m = 800;
c = 0.4;
b = 3600;
umax = 1.6;
umin = -1.4;
acc_comf = 2.5;
gamma = 1.8;
vg = 16;
h = 10;
w = 50;

%% 2.1
vmax = sqrt((b*umax)/(c*(1 + gamma*2)));
acc_max = ((b/m)*umax)/(1 + gamma);
dec_max = ((b/m)*umin)/(1 + gamma*2) - c*vmax^2/m;

%% 2.2
alpha_val = eval(subs(alpha_sol));
beta_val = eval(subs(beta_sol));

idx = alpha_val > 0 & alpha_val < vmax;
alpha_val = alpha_val(idx);
beta_val = beta_val(idx);
%% 2.3

params.g = 9.8;
params.m = 800;
params.c = 0.4;
params.b = 3600;
params.umax = 1.6;
params.umin = -1.4;
params.acc_comf = 2.5;
params.gamma = 1.8;
params.vg = 16;
params.h = 10;
params.w = 50;
params.vmax = vmax;
params.acc_max = acc_max;
params.dec_max = dec_max;
params.alpha = alpha_val;
params.beta = beta_val;

% gear ratio
r = 1;

%%

% Define symbolic variables
syms x(t) v(t) u
theta = 0;

% Define the differential equations
ode1 = diff(x, t) == v;
ode2 = diff(v, t) == ((b/m) * u) / (1 + gamma * r) - (g * sin(theta) * x(t)) - (c/m) * v^2;

% Convert the system of ODEs to MATLAB function handle
odes = [ode2; ode1];
[V, S] = odeToVectorField(odes);

MF = matlabFunction(V, 'vars', {'t', 'Y', 'u'});

tsim = 0:h:50;
v_init = 30;
X0 = [0.1; v_init];
tspan = [tsim(1) tsim(end)];

% Solve the ODEs
[t_og, X_og] = ode45(@(t, Y) MF(t, Y, sin(t)), tspan, X0);
[t_pwa, X_pwa] = ode45(@(t, Y) pwa_friction(t, Y, params, r, theta, sin(t)), tspan, X0);

figure;
sgtitle("Friction Force Comparison");
% Subplot 1: Velocity of follower car
subplot(2, 1, 1);
hold on;
plot(t_og, X_og(:, 2), LineWidth=1.2);
plot(t_og, X_pwa(:, 2), LineWidth=1.2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of the follower car over time');
hold off;
grid on;

% Subplot 2: Position of follower car
subplot(2, 1, 2);
hold on;
plot(t_og, X_og(:, 1), LineWidth=1.2);
plot(t_og, X_pwa(:, 1), LineWidth=1.2);
xlabel('Time (s)');
ylabel('Position (m)');
title('Position of the follower car over time');
hold off;
grid on;
legend("Original", "PWA",'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');


%% 2.4

% Initialize x_init
x_road = 1:5000;
road_height = @(pos) min([((2 * params.h * pos) / params.w); (((params.h * pos) / params.w) + params.h); ...
    (3 * params.h * ones(size(pos))); (((-3 * params.h * pos) / (2 * params.w)) + 9 * params.h)]);

y_road = road_height(x_road);
% Plot the result for visualization
figure;
sgtitle("Road Profile");
subplot(2, 1, 1);
plot(x_road, y_road, LineWidth=1.2);
xlabel('Position');
ylabel('Height (m)');
title('Road Height over Position');
grid on;

subplot(2, 1, 2);
plot(x_road, asin(y_road./x_road), LineWidth=1.2);
xlabel('Position');
ylabel('Slope (rad)');
title('Road Slope over Position');
grid on;

v_init = 40;
x_init = 0.01;
X0 = [x_init; v_init];

tspan = [tsim(1) 50];
[t_model, X_model] = ode45(@(t, Y) pwa_model(t, Y, params, sin(t)), tspan, X0);
theta_prof = asin(road_height(X_model(:, 1)')./X_model(:, 1)');

figure;
sgtitle("Theta Calculation");
% Subplot 1: Velocity of follower car
subplot(3, 1, 1);
hold on;
plot(t_model, X_model(:, 2), LineWidth=1.2);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of the follower car over time');
hold off;
grid on;

% Subplot 2: Position of follower car
subplot(3, 1, 2);
hold on;
plot(t_model, X_model(:, 1), LineWidth=1.2);
xlabel('Time (s)');
ylabel('Position (m)');
title('Position of the follower car over time');
hold off;
grid on;

% subplot(3, 1, 3);
% hold on;
% stairs(t_model, theta_prof, LineWidth=1.2);
% xlabel('Time (s)');
% ylabel('Theta (rad)');
% title('Road Profile');
% hold off;
% grid on;


%% 2.5



%% 2.6 

params.r1 = params.b / (params.m * (1 + params.gamma));
params.r2 = params.b / (params.m * (1 + 2 * params.gamma));
params.Ps = (params.c * params.vmax^2 - params.beta) / (params.vmax - params.alpha);

slopes = [((2 * params.h) / params.w), (((params.h) / params.w)), 0, ((-3 * params.h) / (2 * params.w))];
params.theta1 = (2 * params.h) / params.w;
params.theta2 = params.h / params.w;
params.theta3 = 0;
params.theta4 = (-3 * params.h) / (2 * params.w);
params.Ts = 0.1;

Np = 5; % Prediction Horizon
Nc = 4; % Control Horizon

T_end = 25;
T = 0:params.Ts:T_end;

params.lambda = 0.1;
x_init = 0;
params.xmax = x_init + params.vmax * T_end;
params.xmax = 559.5029;

% Dimensions
dim.nx = 2;
dim.nz = 2;
dim.nd = 5;
dim.nu = 1;
dim.nq = 1;
dim.np = 1;
dim.Np = Np;
dim.Nc = Nc;

% System Dynamics (MLD)
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
    
    f = [0; (-params.Ts*params.g*params.theta4)-(params.Ts / params.m) * ( ...
        params.beta - params.alpha * params.Ps)];

[Ap, Bp, fp] = predmodgen(A1, B1, B2, B3, f, dim);

% Inequality constraints (MLD)
[A_ineq, b_ineq, I_x, I_x0, I_ref] = build_constraints(params, dim);


%% 2.6

x0 = [0.01; 40];
x = zeros(length(T), dim.nx);
x(1,:)  = x0';

for k=1:length(T)-1
    x(k+1,:) = ( MLD_test(params, x(k,:)', sin(k*params.Ts)) )' ;
end

figure;
plot(T, x(:, 2));

%% Simulate MLD uding if conditions 

% x0 = [0.01; 40];
% x = zeros(length(T), dim.nx);
% x(1,:)  = x0';
% 
% for k=1:length(T)-1
%     x(k+1,:) = (simulate_MLD(A1, B1, B2, B3, f, x(k,:)', sin(k*params.Ts), params))' ;
% end
% 
% 
% figure;
% % Subplot 1: Velocity of follower car
% hold on;
% plot(t_model, X_model(:, 2), LineWidth=1.2);
% plot(T, x(:, 2), LineWidth=1.2)
% hold off;
% grid on;

%% 2.7

% vref = 5 * ones(dim.Np, 1);
% 
% x0 = [101; 7];

vref = 23* ones(dim.Np, 1);

x0 = [100; 0.925*params.alpha];

options = optimoptions("intlinprog", "Display", "iter");

obj = [zeros(1, dim.nz*dim.Np), zeros(1, dim.nd*dim.Np), zeros(1, dim.nu*dim.Np), ...
          ones(1, dim.nq*dim.Nc), zeros(1, dim.nq*(dim.Np-dim.Nc)), ...
          params.lambda*ones(1, dim.np*dim.Nc), zeros(1, dim.nq*(dim.Np-dim.Nc))]';

intcon = dim.nz*dim.Np + 1:(dim.nz + dim.nd)*dim.Np;
lb = [-inf*ones(dim.nz*dim.Np,1); zeros(length(intcon),1); -inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];
ub = [inf*ones(dim.nz*dim.Np,1); ones(length(intcon),1); inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];

% System dynamics are used to calculate the inequality constraints
A = (A_ineq + I_x * Bp);
b = (b_ineq - I_x * fp) - (I_x * Ap + I_x0) * x0 + I_ref * vref;
output = intlinprog(obj, intcon, A, b, [], [], lb, ub, [], options);

% z = output(1:dim.nz*dim.Np);
% delta = output(dim.nz*dim.Np + 1: (dim.nz + dim.nd)*dim.Np);
% u = output((dim.nz + dim.nd)*dim.Np + 1: (dim.nz + dim.nd + dim.nu)*dim.Np);
% q = output((dim.nz + dim.nd + dim.nu)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np); 
% p = output((dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq + dim.np)*dim.Np);


%% 2.9

% Reference Velocity
T_ref = [T, T_end:params.Ts:min((T_end+dim.Np)/params.Ts, 30)];

vref = build_vref(T_ref, params);
vref = vref';

x0 = [50.000000000000001; 0.925*params.alpha];

options = optimoptions("intlinprog", "Display", "off");

obj = [zeros(1, dim.nz*dim.Np), zeros(1, dim.nd*dim.Np), zeros(1, dim.nu*dim.Np), ...
          ones(1, dim.nq*dim.Nc), zeros(1, dim.nq*(dim.Np-dim.Nc)), ...
          params.lambda*ones(1, dim.np*dim.Nc), zeros(1, dim.nq*(dim.Np-dim.Nc))]';

intcon = dim.nz*dim.Np + 1:(dim.nz + dim.nd)*dim.Np;
lb = [-inf*ones(dim.nz*dim.Np,1); zeros(length(intcon),1); -inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];
ub = [inf*ones(dim.nz*dim.Np,1); ones(length(intcon),1); inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];

% System dynamics are used to calculate the inequality constraints
A = (A_ineq + I_x * Bp);

u_prev = 0;
x_val = zeros(size(T, 2), dim.nx);
x_val(1, :) = x0';

for k = 1:1:size(T, 2)
% for k = 1
    T(k)

    b = (b_ineq - I_x * fp) - (I_x * Ap + I_x0) * x0 + I_ref * vref(k:k+dim.Np-1);

    [output, fval, flag, msg] = intlinprog(obj, intcon, A, b, [], [], lb, ub, [], options);

    if flag == 1
        z = output(1:dim.nz*dim.Np);
        delta = output(dim.nz*dim.Np + 1: (dim.nz + dim.nd)*dim.Np);
        u = output((dim.nz + dim.nd)*dim.Np + 1: (dim.nz + dim.nd + dim.nu)*dim.Np);
        q = output((dim.nz + dim.nd + dim.nu)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np); 
        p = output((dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq + dim.np)*dim.Np); 
    
        x1_Np = Ap * x0 + Bp * [z; delta; u; q; p] + fp;
        x0 = x1_Np(1:dim.nx, :);

        u_prev = u;
    else
        x1_Np = Ap * x0 + Bp * [z; delta; u_prev; q; p] + fp;
        x0 = x1_Np(1:dim.nx, :);
    end
    x_val(k+1, :) = x0';

end

figure;
plot(T, x_val(1:end-1, 2));


% Cost Function
% J = ones(1,Nc) * q + lambda * ones(1,Nc) * p;

% [X_Np] = predmodgen(X0, u, z, delta, params, Np, Ts); % X_Np = [x;v] (dynamic constraints)
% 
% % Equality Constraints
% X == X_Np;
% delta(:,2) + delta(:,3) + delta(:,4) + delta(:,5) == 1;
% 
% % Inequality Constraints 
% X(:,1) - 50 - (xmax - 50) * (1 - delta(:,2)) <= 0;
% eps(1) - (50 + eps(1)) * delta(:,2) - X(:,1) + 50 <= 0;
% X(:,1) - 100 - (xmax - 100) * (1 - delta(:,3)) <= 0;
% eps(1) - (100 + eps(1)) * delta(:,3) - X(:,1) + 100 <= 0;
% X(:,1) - 200 - (xmax - 200) * (1 - delta(:,4)) <= 0;
% eps(1) - (200 + eps(1)) * delta(:,4) - X(:,1) + 200 <= 0;
% 200 - X(:,1) + eps(1) - 200 * (1 - delta(:,5)) <= 0;
% eps(1) + (200 - xmax - eps(1)) * delta(:,5) + X(:,1) - 200 <= 0;
% 
% X(:,2) - params.vg + eps(1) - (params.vmax - params.vg + eps(1)) * (1 - delta(:,1)) <= 0;
% params.vg - X(:,2) - params.vg * delta(:,1) <= 0;
% X(:,2) - params.alpha - (params.vmax - params.alpha) * (1 - delta(:,6)) <= 0;
% eps(1) - (params.alpha + eps(1)) * delta(:,6) - X(:,2) + params.alpha <= 0;
% 
% z(:,1) - umax * delta(:,1) <= 0;
% umin * delta(:,1) - z(:,1) <= 0;
% z(:,1) - u + umin * (1 - delta(:,1)) <= 0;
% u - umax * (1 - delta(:,1)) - z(:,1) <= 0;
% z(:,2) - params.vmax * delta(:,6) <= 0;
% - z(:,2) <= 0;
% z(:,2) - X(:,2) <= 0;
% X(:,2) - params.vmax * (1 - delta(:,6)) - z(:,2) <= 0;
% 
% u - umax <= 0;
% umin - u <= 0;
% -X(:,2) <= 0;
% X(:,2) - params.vmax <= 0;
% -X(:,1) <= 0;
% X(:,1) - xmax <= 0;
% X(2:end,2) - X(1:end-1,2) - Ts*(params.acc_comf) <= 0;
% -X(2:end,2) + X(1:end-1,2) - Ts*(params.acc_comf) <= 0;
% 
% X(:,2) - vref - q <= 0;
% -X(:,2) + vref - q <= 0;
% u(2:end,2) - u(1:end-1,2) - p <= 0;
% -u(2:end,2) + u(1:end-1,2) - p <= 0;
% u(1) - p <= 0;
% -u(1) - p <= 0;
% 


%% Function 

function dx = pwa_friction(t, x, params, r, theta, u)
    dx = zeros(2, 1);
    pos = x(1);
    vel = x(2);

    p1 = (params.beta/params.alpha)*vel;
    p2 = ((params.c*params.vmax^2-params.beta) / (params.vmax-params.alpha))*(vel - params.alpha) + params.beta;

    if vel <= params.alpha
        acc = ((params.b/params.m) * u) / (1 + params.gamma * r) - (params.g * sin(theta)) - (p1 / params.m);
    elseif vel > params.alpha
        acc = ((params.b/params.m) * u) / (1 + params.gamma * r) - (params.g * sin(theta)) - (p2 / params.m);
    else
        disp("Error");
    end

    dx = [vel; acc];
end


function [dx] = pwa_model(t, x, params, u)
    dx = zeros(2, 1);
    pos = x(1);
    vel = x(2);

    slopes = [((2 * params.h) / params.w), (((params.h) / params.w)), 0, ((-3 * params.h) / (2 * params.w))];
    [y, idx] = min([((2 * params.h * pos) / params.w), (((params.h * pos) / params.w) + params.h), (3 * params.h), (((-3 * params.h * pos) / (2 * params.w)) + 9 * params.h)]);
    theta = slopes(idx);
    
    p1 = (params.beta / params.alpha) * vel;
    p2 = ((params.c * params.vmax^2 - params.beta) / (params.vmax - params.alpha)) * (vel - params.alpha) + params.beta;

    if vel >= 0 && vel < params.vg
        r = 1;
        acc = ((params.b / params.m) * u) / (1 + params.gamma * r) - (params.g * (theta)) - (p1 / params.m);
    elseif vel >= params.vg && vel <= params.alpha
        r = 2;
        acc = ((params.b / params.m) * u) / (1 + params.gamma * r) - (params.g * (theta)) - (p1 / params.m);
    elseif vel > params.alpha
        r = 2;
        acc = ((params.b / params.m) * u) / (1 + params.gamma * r) - (params.g * (theta)) - (p2 / params.m);
    else
        disp('Error');
        return;
    end

    % acc = max(min(acc, params.acc_max), params.dec_max);
    % vel = max(min(vel, params.vmax), 0)

    dx = [vel; acc];
end

% function [X_Np] = predmodgen(X0, u, z, delta, params, Np, Ts)
%     % X0 (initial states) = 1x2 (x0, v0)
%     % size(u) = Npx1
%     % size(z) = Npx2
%     % size(delta) = Npx6
%     % size(X_Np) = Npx2
% 
%     X_Np = zeros(Np,2);
%     X_Np(1,:) = X0;
% 
%     for i = 1:Np
%         X_Np(i+1,1) = X_Np(i,1) + Ts * X_Np(i,2);
%         X_Np(i+1,2) = X_Np(i,2) + Ts* ((params.r1 - params.r2) * z(i,1) ...
%                     + params.r2 * u - (1 / params.m) * ((z(i,2) * (params.beta / params.alpha) - params.Ps) ...
%                     + params.Ps * (X_Np(i,2) - params.alpha) + params.beta + delta(i,6) * (params.Ps * params.alpha - params.beta)) ...
%                     - params.g * (delta(i,2) * params.theta1 + delta(i,3) * params.theta2 + delta(i,5) * params.theta4));
%     end
% 
% end


function [Ap, Bp, fp] = predmodgen(A1, B1, B2, B3, f, dim)

    % Prediction matrices generation
    % This function computes the prediction matrices to be used in the
    % optimization problem

    % Prediction matrix from initial state
    % This will give vector x(1) to x(Np)
    Ap = zeros(dim.nx * dim.Np, dim.nx);
    for k = 0:dim.Np-1
        Ap(k*dim.nx+1:(k+1)*dim.nx, :) = A1^(k+1);
    end
    
    B_1p = zeros(dim.nx * dim.Np, dim.nz * dim.Np); % Prediction matrix from zp
    B_2p = zeros(dim.nx * dim.Np, dim.nd * dim.Np); % Prediction matrix from delta_p
    B_3p = zeros(dim.nx * dim.Np, dim.nu * dim.Np); % Prediction matrix from up
    fp = zeros(dim.nx * dim.Np, 1); % Constant
    fp(1:dim.nx, 1) = f;

    for k = 0:dim.Np-1
        for i = 0:k
            B_1p(k*dim.nx+1:(k+1)*dim.nx, i*dim.nz+1:(i+1)*dim.nz) = A1^(k-i) * B1;

            B_2p(k*dim.nx+1:(k+1)*dim.nx, i*dim.nd+1:(i+1)*dim.nd) = A1^(k-i) * B2;

            B_3p(k*dim.nx+1:(k+1)*dim.nx, i*dim.nu+1:(i+1)*dim.nu) = A1^(k-i) * B3;
        end

        fp((k+1)*dim.nx+1:(k+2)*dim.nx, 1) = A1 * (fp(k*dim.nx+1:(k+1)*dim.nx, 1)) + f;
    end

    fp = fp(1:end-dim.nx, 1);

    Bp = [B_1p, B_2p, B_3p, zeros(dim.nx * dim.Np, dim.Np), zeros(dim.nx * dim.Np, dim.Np)];

end

function vref = build_vref(t, params)
    % Initialize the output
    vref = zeros(size(t));

    % Define the intervals
    interval1 = (t >= 0 & t <= 3);
    interval2 = (t > 3 & t <= 9);
    interval3 = (t > 9 & t <= 15);
    interval4 = (t > 15 & t <= 18);
    interval5 = (t > 18 & t <= 21);
    interval6 = (t > 21 & t <= 30);

    % Apply the function for each interval
    vref(interval1) = 0.85 * params.alpha;
    vref(interval2) = 1.2 * params.alpha;
    vref(interval3) = 1.2 * params.alpha - (1/12) * params.alpha * (t(interval3) - 9);
    vref(interval4) = 0.7 * params.alpha;
    vref(interval5) = 0.7 * params.alpha + (4/15) * params.alpha * (t(interval5) - 18);
    vref(interval6) = 0.9 * params.alpha;
end



function x_next = simulate_MLD(A1, B1, B2, B3, f, x0, u, params)
    pos = x0(1);
    vel = x0(2);

    delta_1 = 0; delta_2 = 0; delta_3 = 0; delta_4 = 0; delta_5 = 0;
    
    if (pos >= 0) && (pos <= 50)
        delta_2 = 1; delta_3 = 1; delta_4 = 1; delta_5 = 1;
    elseif (pos >= 50 + eps(1)) && (pos <= 100)
        delta_3 = 1;
    elseif (pos >= 100 + eps(1)) && (pos <= 200)
        delta_4 = 1;
    elseif (pos >= 200 + eps(1))
        delta_5 = 1;
    else
        disp("Error");
    end

    if (vel >= 0) && (vel <= params.vg - eps(1))
        delta_1 = 1;
    elseif (vel >= params.vg) && (vel <= params.alpha)
        delta_6 = 1;
    elseif (vel > params.alpha) && (vel <= params.vmax)
        disp("vel btw alpha and vmax")
    elseif (vel > params.vmax)
        disp("vel greater than vmax")
    end

    z_1 = delta_1 * u;
    z_2 = delta_6 * vel;

    delta = [delta_1; delta_2; delta_3; delta_4; delta_5; delta_6];
    z = [z_1; z_2];

    x_next = A1 * x0 + B1 * z + B2 * delta + B3 * u + f;
end