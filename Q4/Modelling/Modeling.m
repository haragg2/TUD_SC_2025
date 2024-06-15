clear all;
close all;
clc;

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%%

% Parameters
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

params.Ts = 0.1;

%For Q2.3
% gear ratio
r = 1;

syms v vmax alpha beta

p1 = (beta/alpha)*v;
p2 = ((params.c*vmax^2-beta) / (vmax-alpha))*(v - alpha) + beta;

a1_int = @(v) (params.c * v^2 - p1)^2;
a1 =  simplify(int(a1_int, v, 0, alpha));

a2_int = @(v) (params.c * v^2 - p2)^2;
a2 = simplify(int(a2_int, v, alpha, vmax));

pa = a1+a2;
dp_alpha = diff(pa, alpha);
dp_beta = diff(pa, beta);

% Solve the system of equations
solutions = solve([dp_alpha == 0, dp_beta == 0], [alpha, beta]);

alpha_sol = solutions.alpha;
beta_sol = solutions.beta;



%% 2.1
vmax = sqrt((params.b * params.umax)/(params.c * (1 + params.gamma * 2)));
acc_max = ((params.b / params.m) * params.umax)/(1 + params.gamma);
dec_max = ((params.b / params.m) * params.umin)/(1 + params.gamma * 2) - params.c * vmax^2/ params.m;

params.vmax = vmax;
params.acc_max = acc_max;
params.dec_max = dec_max;

%% 2.2
alpha_val = eval(subs(alpha_sol));
beta_val = eval(subs(beta_sol));

idx = alpha_val > 0 & alpha_val < vmax;
alpha_val = alpha_val(idx);
beta_val = beta_val(idx);

params.alpha = alpha_val;
params.beta = beta_val;

% Define piecewise function for P(v)
P = piecewise(v <= params.alpha, (params.beta / params.alpha) * v, ...
              params.alpha < v & v <= params.vmax, ...
              (params.c*params.vmax^2 - params.beta) / (params.vmax - params.alpha) * (v - params.alpha) + params.beta);

% Generate values of v from 0 to vmax
v_values = linspace(0, params.vmax, 1000);

% Evaluate P(v) for the generated v values
P_pwa = double(subs(P, v, v_values));
P_nl = params.c * v_values.^2;

% Generate values of v from 0 to vmax
v = linspace(0, params.vmax, 1000);

% Plotting
figure;
hold on;
plot(v_values, P_nl, 'LineWidth', 1.2);
plot(v_values, P_pwa, 'LineWidth', 1.2);
plot([params.alpha, params.alpha], [0, params.beta], 'k--');
plot([0, params.alpha], [params.beta, params.beta], 'k--');
text(params.alpha, -30, '$$\alpha$$', 'HorizontalAlignment', 'center');
text(-1, params.beta + 30, '$$\beta$$', 'HorizontalAlignment', 'center');
hold off;
xlabel('Velocity ($m/s$)');
ylabel('$F$');
legend("Non-linear (V)", "PWA (P)");
title('Friction Curve');
grid on;

%% 2.3

% Define symbolic variables
syms x(t) v(t) u
theta = 0;

% Define the differential equations
ode1 = diff(x, t) == v;
ode2 = diff(v, t) == ((params.b/params.m) * u) / (1 + params.gamma * r) - (params.g * sin(theta) * x(t)) - (params.c / params.m) * v^2;

% Convert the system of ODEs to MATLAB function handle
odes = [ode2; ode1];
[V, S] = odeToVectorField(odes);

MF = matlabFunction(V, 'vars', {'t', 'Y', 'u'});

T_end = 50;
T = 0:params.Ts:T_end;
v_init = 30;
x0 = [0.1; v_init];
tspan = [T(1) T(end)];

% Solve the ODEs
[t_og, X_og] = ode45(@(t, Y) MF(t, Y, sin(t)), tspan, x0);
[t_pwa, X_pwa] = ode45(@(t, Y) pwa_friction(t, Y, params, r, theta, sin(t)), tspan, x0);

figure;
sgtitle("Friction Force Comparison");
subplot(2, 1, 1);
hold on;
plot(t_og, X_og(:, 1), LineWidth=1.2);
plot(t_og, X_pwa(:, 1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Position ($m$)');
title('Position over time');
hold off;
grid on;

subplot(2, 1, 2);
hold on;
plot(t_og, X_og(:, 2), LineWidth=1.2);
plot(t_og, X_pwa(:, 2), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
title('Velocity over time');
hold off;
grid on;
legend("NL Model", "PWA Model",'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');


%% 2.4

% Initialize x_init
x_road = 1:250;
slopes = [((2 * params.h) / params.w); (params.h / params.w); 0; ((-3 * params.h) / (2 * params.w))];

road_height = @(pos) min([((2 * params.h * pos) / params.w); (((params.h * pos) / params.w) + params.h); ...
    (3 * params.h * ones(size(pos))); (((-3 * params.h * pos) / (2 * params.w)) + 9 * params.h)]);

y_road = road_height(x_road);

% Plot the result for visualization
figure;
sgtitle("Road Profile");
subplot(2, 1, 1);
plot(x_road, y_road, LineWidth=1.2);
xlabel('Position ($m$)');
ylabel('Height ($m$)');
title('Road Height over Position');
grid on;

subplot(2, 1, 2);
plot(x_road(2:end), asin((y_road(2:end) - y_road(1:end-1)) ./ (x_road(2:end) - x_road(1:end-1))), LineWidth=1.2);
xlabel('Position ($m$)');
ylabel('Slope (rad)');
title('Road Slope over Position');
grid on;

v_init = 40;
x_init = 10;
x0 = [x_init; v_init];
T_end = 40;

tspan = [T(1) T_end];
[t_pwa, x_pwa_cont] = ode45(@(t, Y) pwa_model(t, Y, params, 0.8 * cos(t)), tspan, x0);
[t_nl, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, 0.8 * cos(t), params), tspan, x0);
theta_prof = asin(road_height(x_pwa_cont(:, 1)')./x_pwa_cont(:, 1)');

figure;
sgtitle("PWA Model Simulation");
subplot(2, 1, 1);
hold on;
plot(t_nl, x_nl(:, 1), LineWidth=1.2);
plot(t_pwa, x_pwa_cont(:, 1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Position ($m$)');
title('Position over time');
hold off;
grid on;

subplot(2, 1, 2);
hold on;
plot(t_nl, x_nl(:, 2), LineWidth=1.2);
plot(t_pwa, x_pwa_cont(:, 2), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
title('Velocity over time');
hold off;
grid on;
legend("NL Model", "PWA Model",'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');

figure;
theta_pwa = zeros(length(x_pwa_cont(:, 1)), 1);
for k=1:1:length(theta_pwa)
    pos = x_pwa_cont(k, 1);
    [y, idx] = min([((2 * params.h * pos) / params.w), (((params.h * pos) / params.w) + params.h), (3 * params.h), (((-3 * params.h * pos) / (2 * params.w)) + 9 * params.h)]);
    theta_pwa(k) = slopes(idx);
end
plot(t_pwa, theta_pwa, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Theta (rad)');
title('Road Profile ($\theta$)');
grid on;

%% 2.5

T_end = 40;
T = 0:params.Ts:T_end;

% Initial condition same as 2.4
x_pwa_dis = zeros(length(T), 2);
x_pwa_dis(1, :) = x0;

for k=1:1:length(T)-1
    u = 0.8 * cos(k * params.Ts);
    dx = pwa_model(k, x_pwa_dis(k, :), params, u);

    x_pwa_dis(k+1, :) = (x_pwa_dis(k, :)' + dx * params.Ts)';
end

figure;
sgtitle("Forward Euler Discretization of PWA model");
subplot(2, 1, 1);
hold on;
plot(t_pwa, x_pwa_cont(:, 1), LineWidth=1.2);
plot(T, x_pwa_dis(:, 1), LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Position ($m$)');
grid on;

subplot(2, 1, 2);
hold on;
plot(t_pwa, x_pwa_cont(:, 2), LineWidth=1.2);
plot(T, x_pwa_dis(:, 2), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
hold off;
grid on;
legend("PWA Continuous", "PWA Discrete",'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');


%% 2.6

params.r1 = params.b / (params.m * (1 + params.gamma));
params.r2 = params.b / (params.m * (1 + 2 * params.gamma));
params.Ps = (params.c * params.vmax^2 - params.beta) / (params.vmax - params.alpha);

slopes = [((2 * params.h) / params.w), (((params.h) / params.w)), 0, ((-3 * params.h) / (2 * params.w))];
params.theta1 = (2 * params.h) / params.w;
params.theta2 = params.h / params.w;
params.theta3 = 0;
params.theta4 = (-3 * params.h) / (2 * params.w);

params.lambda = 0.001;
params.xmax = 0 + params.vmax * T_end;

% Dimensions
dim.nx = 2;
dim.nz = 2;
dim.nd = 5;
dim.nu = 1;
dim.nq = 1;
dim.np = 1;
dim.Np = 5;
dim.Nc = 4;

% To verify the MLD model
x0 = [0.01; 40];
x = zeros(length(T), dim.nx);
x(1,:)  = x0';

for k=1:length(T)-1
    x(k+1,:) = ( MLD_test(params, x(k,:)', sin(k*params.Ts)) )' ;
end

figure;
sgtitle("MLD vs PWA Model Simulation");
subplot(2, 1, 1);
hold on;
plot(T, x_pwa_dis(:, 1), LineWidth=1.2);
plot(T, x(:, 1), LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Position ($m$)');
grid on;

subplot(2, 1, 2);
hold on;
plot(T, x_pwa_dis(:, 2), LineWidth=1.2);
plot(T, x(:, 2), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
hold off;
grid on;
legend("PWA Discrete", "MLD Model",'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');

%% 2.7

% Update Prediction and Control Horizon
dim.Np = 4;
dim.Nc = 2;

% Vref for the follower car
vref = 40.2 * ones(dim.Np, 1);
x0 = [100; 40];

% Update max distance
params.xmax = x0(1) + params.vmax * T_end;

% Update lambda for the objective computation
params.lambda = 0.001;

% Get the objective function and prediction matrices
[pred, obj] = build_obj(params, dim);

% System dynamics are used to calculate the inequality constraints
A = (obj.A_ineq + obj.I_x * pred.Bp);
b = (obj.b_ineq - obj.I_x * pred.fp) - (obj.I_x * pred.Ap + obj.I_x0) * x0 + obj.I_ref * vref;

% Solve the MLD MPC for an arbitrary k
% output = intlinprog(obj.func, obj.intcon, A, b, obj.Aeq, obj.beq, obj.lb, obj.ub, [], obj.options);

% GLPK
ctype = [repmat('U', 1, size(A, 1), 1) repmat('S', 1, size(obj.Aeq, 1), 1)];
vartype = [repmat('C', 1, dim.nz*dim.Np) repmat('I', 1, dim.nd*dim.Np) repmat('C', 1, dim.nu*dim.Np) ...
            repmat('C', 1, dim.nq*dim.Np) repmat('C', 1, dim.np*dim.Np)];

[output_glpk, fval, flag, msg] = glpk(obj.func, [A; obj.Aeq], [b; obj.beq], obj.lb, obj.ub, ctype, vartype, 1);

optimal_control = output_glpk((dim.nz + dim.nd)*dim.Np + 1: (dim.nz + dim.nd + dim.nu)*dim.Np);

% Simulate the MLD system with one control sequence 
x_val_Np = pred.Ap * x0 + pred.Bp * output_glpk + pred.fp;
x_val_Np = reshape(x_val_Np,  [2, dim.Np])';


figure(7);
sgtitle("MPC Input Sequence for arbitrary step $k$");
subplot(2, 1, 1);
hold on;
plot(T(1:dim.Np), optimal_control, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
legendu = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legendu1 = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
legend(legendu);
% legend(legendu, legendu1);
grid on;
hold off;

subplot(2, 1, 2);
hold on;
plot(T(1:dim.Np), x_val_Np(:, 2), LineWidth=1.2);
% plot(T(1:dim.Np), vref(1:dim.Np), '--', LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
legendv = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legendv1 = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
legend(legendv,'Reference Velocity');
% legend(legendv,'Reference Velocity', legendv1);
grid on;

%% 2.8

dim.Np = 4;
dim.Nc = 2;

% Simulation time
T_end = 5;
T = 0:params.Ts:T_end;

% Initial condition
x0 = [100; 40];

% Update max distance
params.xmax = x0(1) + params.vmax * T_end;

% Reference Velocity
T_ref = [T, T_end:params.Ts:(T_end+dim.Np)/params.Ts];
vref = 40.2*ones(length(T_ref), 1);

% Update lambda for the objective computation
params.lambda = 0.001;

% Get the objective function and prediction matrices
[pred, obj] = build_obj(params, dim);

% System dynamics are used to calculate the inequality constraints
A = (obj.A_ineq + obj.I_x * pred.Bp);

u_prev = 0;
x_val = zeros(size(T, 2), dim.nx);
x_val(1, :) = x0';

optimal_control = zeros(length(T)-1, 1);

ctype = [repmat('U', 1, size(A, 1), 1) repmat('S', 1, size(obj.Aeq, 1), 1)];
vartype = [repmat('C', 1, dim.nz*dim.Np) repmat('I', 1, dim.nd*dim.Np) repmat('C', 1, dim.nu*dim.Np) ...
            repmat('C', 1, dim.nq*dim.Np) repmat('C', 1, dim.np*dim.Np)];

for k = 1:1:size(T, 2)-1
    tspan = [T(k) T(k+1)];

    % This inequality depends on the initial x0
    b = (obj.b_ineq - obj.I_x * pred.fp) - (obj.I_x * pred.Ap + obj.I_x0) * x_val(k, :)' + obj.I_ref * vref(k:k+dim.Np-1);

    % [output, fval, flag, msg] = intlinprog(obj.func, obj.intcon, A, b, obj.Aeq, obj.beq, obj.lb, obj.ub, [], obj.options);
    [output, fval, flag, msg] = glpk(obj.func, [A; obj.Aeq], [b; obj.beq], obj.lb, obj.ub, ctype, vartype, 1);

    % if flag == 1
    if flag == 5
        z = output(1:dim.nz*dim.Np);
        delta = output(dim.nz*dim.Np + 1: (dim.nz + dim.nd)*dim.Np);
        u = output((dim.nz + dim.nd)*dim.Np + 1: (dim.nz + dim.nd + dim.nu)*dim.Np);
        q = output((dim.nz + dim.nd + dim.nu)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np); 
        p = output((dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq + dim.np)*dim.Np); 

        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u(1), params), tspan, x0);
        x0 = x_nl(end, :)';

        u_prev = u(1);
    else
        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u_prev, params), tspan, x0);
        x0 = x_nl(end, :)';
    end

    optimal_control(k, :) = u_prev;
    x_val(k+1, :) = x0;

end

% Figures
figure(8);
sgtitle("MPC Closed-loop Response");
subplot(2, 1, 1);
hold on;
plot(T(1:end-1), optimal_control, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
% legendu = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legend(legendu);
legendu1 = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
legend(legendu, legendu1);
grid on;
hold off;

subplot(2, 1, 2);
hold on;
plot(T, x_val(:, 2), LineWidth=1.2);
% plot(T, vref(1:length(T)), '--', LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
% legendv = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legend(legendv,'Reference Velocity');
legendv1 = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
legend(legendv,'Reference Velocity', legendv1);
grid on;

%% 2.9

% Simulation time
T_end = 25;
T = 0:params.Ts:T_end;

% Reference Velocity
T_ref = [T, T_end:params.Ts:min((T_end+dim.Np)/params.Ts, 30)];
vref = build_vref(T_ref, params);
vref = vref';

% Initial condition
x0 = [101; 0.925*params.alpha];

% Update Prediction and Control Horizon
dim.Np = 5;
dim.Nc = 4;

% Update max distance using initial condition and vmax
params.xmax = x0(1) + params.vmax * T_end;

% Update lambda for the objective computation
params.lambda = 0.1;

% Get the objective function and prediction matrices
[pred, obj] = build_obj(params, dim);

% System dynamics are used to calculate the inequality constraints
A = (obj.A_ineq + obj.I_x * pred.Bp);

u_prev = 0;
x_val = zeros(size(T, 2), dim.nx);
x_val(1, :) = x0';

optimal_control = zeros(length(T)-1, 1);

ctype = [repmat('U', 1, size(A, 1), 1) repmat('S', 1, size(obj.Aeq, 1), 1)];
vartype = [repmat('C', 1, dim.nz*dim.Np) repmat('I', 1, dim.nd*dim.Np) repmat('C', 1, dim.nu*dim.Np) ...
            repmat('C', 1, dim.nq*dim.Np) repmat('C', 1, dim.np*dim.Np)];

% To measure the computation time
total_time = 0;
num_iterations = size(T, 2) - 1;
iteration_times_imp = zeros(1, num_iterations);

for k = 1:1:num_iterations
    % Start timing
    tic;

    tspan = [T(k) T(k+1)];

    % This inequality depends on the initial x0
    b = (obj.b_ineq - obj.I_x * pred.fp) - (obj.I_x * pred.Ap + obj.I_x0) * x0 + obj.I_ref * vref(k:k+dim.Np-1);

    % [output, fval, flag, msg] = intlinprog(obj.func, obj.intcon, A, b, obj.Aeq, obj.beq, obj.lb, obj.ub, [], obj.options);
    [output, fval, flag, msg] = glpk(obj.func, [A; obj.Aeq], [b; obj.beq], obj.lb, obj.ub, ctype, vartype, 1);

    iteration_time = toc; % Stop timing
    iteration_times_imp(k) = iteration_time; % Store the computation time for this iteration
    total_time = total_time + iteration_time; % Accumulate the total time

    % if flag == 1
    if flag == 5
        z = output(1:dim.nz*dim.Np);
        delta = output(dim.nz*dim.Np + 1: (dim.nz + dim.nd)*dim.Np);
        u = output((dim.nz + dim.nd)*dim.Np + 1: (dim.nz + dim.nd + dim.nu)*dim.Np);
        q = output((dim.nz + dim.nd + dim.nu)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np); 
        p = output((dim.nz + dim.nd + dim.nu + dim.nq)*dim.Np + 1: (dim.nz + dim.nd + dim.nu + dim.nq + dim.np)*dim.Np); 

        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u(1), params), tspan, x0);
        x0 = x_nl(end, :)';

        u_prev = u(1);
    else
        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u_prev, params), tspan, x0);
        x0 = x_nl(end, :)';
    end

    optimal_control(k, :) = u_prev;
    x_val(k+1, :) = x0;
end

avg_time_per_itr = total_time / num_iterations;
fprintf(['Average time to compute control input with Implicit MPC: %.4f seconds' ...
    ' with Np = %d, Nc = %d\n'], avg_time_per_itr, dim.Np, dim.Nc);

% Figures
figure;
gTitle = sprintf('Closed-loop States Evolution with Implicit MPC\n$N_p=%d$, $N_c=%d$, $\\lambda=%.1f$', dim.Np, dim.Nc, params.lambda);
sgtitle(gTitle);
subplot(3, 1, 1);
plot(x_val(:, 1), x_val(:, 2), LineWidth=1.2);
xlabel('Position ($m$)');
ylabel('Velocity ($m/s$)');
title("State Evolution");
grid on;
subplot(3, 1, 2);
plot(T, x_val(:, 1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Position ($m$)');
title("Position over Time");
grid on;
subplot(3, 1, 3);
hold on;
plot(T, x_val(:, 2), LineWidth=1.2);
plot(T, vref(1:length(T)), '--', LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
legend("$v$", "$v_{ref}$");
title("Velocity over Time");
grid on;

figure;
gTitle = sprintf('Closed-loop Response with Implicit MPC\n$N_p=%d$, $N_c=%d$, $\\lambda=%.1f$', dim.Np, dim.Nc, params.lambda);
sgtitle(gTitle);
subplot(4, 1, 1);
plot(T(1:end-1), optimal_control, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
title("Throttle Input over Time");
grid on;
subplot(4, 1, 2);
plot(T(2:end-1), optimal_control(2:end) - optimal_control(1:end-1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
title("Throttle Input Change over Time");
grid on;
subplot(4, 1, 3);
plot(T, x_val(:, 2) - vref(1:length(T)), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
title("Velocity Error over Time");
grid on;
subplot(4, 1, 4);
plot(T(2:end), (x_val(2:end, 2) - x_val(1:end-1, 2))/params.Ts, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Acceleration ($m/s^2$)');
title("Acceleration over Time");
grid on;

% Plot the computation time vs iteration
figure;
plot(T(1:end-1), iteration_times_imp, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Computation Time ($s$)');
title(sprintf('Computation Time for Control Input with Implicit MPC\n$N_p=%d$, $N_c=%d$', dim.Np, dim.Nc));
grid on;

% figure;
% sgtitle("Implicit MPC Performance Comparison");
% subplot(2, 1, 1);
% hold on;
% plot(T_54(1:end-1), u_54, LineWidth=1.2);
% plot(T(1:end-1), optimal_control, LineWidth=1.2);
% hold off;
% xlabel('Time ($s$)');
% ylabel('Throttle Input');
% title("Throttle Input over Time");
% leg_v1 = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', 5, 4, params.lambda);
% leg_v2 = sprintf('$u$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legend(leg_v1, leg_v2);
% grid on;
% subplot(2, 1, 2);
% hold on;
% plot(T_54, vel_54, LineWidth=1.2);
% plot(T, x_val(:, 2), LineWidth=1.2);
% plot(T, vref(1:length(T)), '--', LineWidth=1.2);
% hold off;
% xlabel('Time ($s$)');
% ylabel('Velocity ($m/s$)');
% leg_v1 = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', 5, 4, params.lambda);
% leg_v2 = sprintf('$v$ at $N_p=%d$, $N_c=%d$, $\\lambda=%.3f$', dim.Np, dim.Nc, params.lambda);
% legend(leg_v1, leg_v2, "$v_{ref}$");
% title("Velocity over Time");
% grid on;

%% 2.10

dim.Np = 2;
dim.Nc = 2;

x0_var = sdpvar(dim.nx, 1);
z_var = sdpvar(dim.nz * dim.Np, 1);
d_var = binvar(dim.nd * dim.Np, 1);
u_var = sdpvar(dim.nu * dim.Np, 1);
q_var = sdpvar(dim.nq * dim.Np, 1);
p_var = sdpvar(dim.np * dim.Np, 1);
vref_var = sdpvar(dim.Np, 1);

% Get the objective function and prediction matrices
[pred, obj] = build_obj(params, dim);

% System dynamics are used to calculate the inequality constraints
A = (obj.A_ineq + obj.I_x * pred.Bp);
b = (obj.b_ineq - obj.I_x * pred.fp) - (obj.I_x * pred.Ap + obj.I_x0) * x0_var + obj.I_ref * vref_var;

X = [z_var; d_var; u_var; q_var; p_var];
objective = obj.func' * X;
constraints = [A * X <= b; obj.Aeq * X == obj.beq];

tic;
[sol, diagn, Z, Valuefcn, Optimizer] = solvemp(constraints, objective, [], [x0_var; vref_var], u_var(1));
explicit_mpc_comp_time = toc;

fprintf(['Total time to compute Explicit MPC: %.4f seconds' ...
    ' with Np = %d, Nc = %d\n'], explicit_mpc_comp_time, dim.Np, dim.Nc);

% Initial condition
x0_ex = [101; 0.925*params.alpha];

u_prev = 0;
x_val_ex = zeros(size(T, 2), dim.nx);
x_val_ex(1, :) = x0_ex';

optimal_control_ex = zeros(length(T)-1, 1);

% To measure the computation time
total_time = 0;
num_iterations = size(T, 2) - 1;
iteration_times_exp = zeros(1, num_iterations);

for k = 1:1:num_iterations
    % Start timing
    tic;

    tspan = [T(k) T(k+1)];

    assign(x0_var, x0_ex);
    assign(vref_var, vref(k:k+dim.Np-1));
    u = value(Optimizer);

    iteration_time = toc; % Stop timing
    iteration_times_exp(k) = iteration_time; % Store the computation time for this iteration
    total_time = total_time + iteration_time; % Accumulate the total time

    if isnan(u)
        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u_prev, params), tspan, x0_ex);
        x0_ex = x_nl(end, :)';
    else
        [~, x_nl] = ode45(@(t, Y) NL_Dynamics(t, Y, u, params), tspan, x0_ex);
        x0_ex = x_nl(end, :)';
        u_prev = u;
    end

    optimal_control_ex(k, :) = u_prev;
    x_val_ex(k+1, :) = x0_ex;
end

avg_time_per_itr = total_time / num_iterations;
fprintf(['Average time to compute control input with Explicit MPC: %.4f seconds' ...
    ' with Np = %d, Nc = %d\n'], avg_time_per_itr, dim.Np, dim.Nc);

% Figures
figure;
gTitle = sprintf('Closed-loop States Evolution with Explicit MPC\n$N_p=%d$, $N_c=%d$, $\\lambda=%.1f$', dim.Np, dim.Nc, params.lambda);
sgtitle(gTitle);%, 'Interpreter', 'latex');

subplot(3, 1, 1);
plot(x_val_ex(:, 1), x_val_ex(:, 2), LineWidth=1.2);
xlabel('Position ($m$)');
ylabel('Velocity ($m/s$)');
title("State Evolution");
grid on;
subplot(3, 1, 2);
plot(T, x_val_ex(:, 1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Position ($m$)');
title("Position over Time");
grid on;
subplot(3, 1, 3);
hold on;
plot(T, x_val_ex(:, 2), LineWidth=1.2);
plot(T, vref(1:length(T)), '--', LineWidth=1.2);
hold off;
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
legend("$v$", "$v_{ref}$");
title("Velocity over Time");
grid on;

figure;
gTitle = sprintf('Closed-loop Response with Explicit MPC\n$N_p=%d$, $N_c=%d$, $\\lambda=%.1f$', dim.Np, dim.Nc, params.lambda);
sgtitle(gTitle);
subplot(4, 1, 1);
plot(T(1:end-1), optimal_control_ex, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
title("Throttle Input over Time");
grid on;
subplot(4, 1, 2);
plot(T(2:end-1), optimal_control_ex(2:end) - optimal_control_ex(1:end-1), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Throttle Input');
title("Throttle Input Change over Time");
grid on;
subplot(4, 1, 3);
plot(T, x_val_ex(:, 2) - vref(1:length(T)), LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Velocity ($m/s$)');
title("Velocity Error over Time");
grid on;
subplot(4, 1, 4);
plot(T(2:end), (x_val_ex(2:end, 2) - x_val_ex(1:end-1, 2))/params.Ts, LineWidth=1.2);
xlabel('Time ($s$)');
ylabel('Acceleration ($m/s^2$)');
title("Acceleration over Time");
grid on;

%% Function

function dx = NL_Dynamics(t, x, u, params)
    pos = x(1); vel = x(2);

    x_dot = vel;

    slopes = [((2 * params.h) / params.w), (((params.h) / params.w)), 0, ((-3 * params.h) / (2 * params.w))];
    [y, idx] = min([((2 * params.h * pos) / params.w), (((params.h * pos) / params.w) + params.h), (3 * params.h), (((-3 * params.h * pos) / (2 * params.w)) + 9 * params.h)]);
    theta = slopes(idx);

    if vel >= 0 && vel < params.vg
        r = 1;
    else
        r = 2;
    end

    v_dot = ((params.b / params.m) * u) / (1 + params.gamma * r) - (params.g * sin(theta)) - (params.c / params.m) * vel^2;

    dx = [x_dot; v_dot];
end


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
    end

    dx = [vel; acc];
end


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
    end

    if (vel >= 0) && (vel <= params.vg - eps(1))
        delta_1 = 1;
    elseif (vel >= params.vg) && (vel <= params.alpha)
        delta_6 = 1;
    elseif (vel > params.alpha) && (vel <= params.vmax)
        disp("vel between alpha and vmax")
    elseif (vel > params.vmax)
        disp("vel greater than vmax")
    end

    z_1 = delta_1 * u;
    z_2 = delta_6 * vel;

    delta = [delta_1; delta_2; delta_3; delta_4; delta_5; delta_6];
    z = [z_1; z_2];

    x_next = A1 * x0 + B1 * z + B2 * delta + B3 * u + f;
end


function [pred, obj] = build_obj(params, dim)
    
    % System Dynamics (MLD)
    pred.A1 = [1, params.Ts; 0, (1 - (params.Ts * params.Ps) / params.m)];
    
    pred.B1 = [0, 0; 
        params.Ts * (params.r1 - params.r2), -(params.Ts / ...
        params.m) * ((params.beta/params.alpha) - params.Ps)];
    
    pred.B2 = [zeros(1, 5);
        0, -(params.Ts * params.g * params.theta2), -( ...
        params.Ts * params.g * params.theta2), ( ...
        params.Ts * params.g * params.theta4), ...
        -(params.Ts / params.m) * (params.alpha * params.Ps - params.beta)];
    
    pred.B3 = [0; params.Ts * params.r2];
    
    pred.f = [0; (-params.Ts*params.g*params.theta4)-(params.Ts / params.m) * ( ...
        params.beta - params.alpha * params.Ps)];
    
    % Get prediction matrices
    [pred.Ap, pred.Bp, pred.fp] = predmodgen(pred.A1, pred.B1, pred.B2, pred.B3, pred.f, dim);
    
    % Inequality constraints (MLD)
    [obj.A_ineq, obj.b_ineq, obj.I_x, obj.I_x0, obj.I_ref] = build_constraints(params, dim);
    
    % Equality constraints
    obj.Aeq = [zeros(dim.Np-dim.Nc, dim.nz*dim.Np), zeros(dim.Np-dim.Nc, dim.nd*dim.Np), ...
           zeros(dim.Np-dim.Nc, dim.Nc-1), -1*ones(dim.Np-dim.Nc, 1), eye(dim.Np-dim.Nc),...
           zeros(dim.Np-dim.Nc, dim.nq*dim.Np), zeros(dim.Np-dim.Nc, dim.np*dim.Np)];
    
    obj.beq = zeros(dim.Np-dim.Nc, 1);
    
    obj.options = optimoptions("intlinprog", "Display", "off");
    obj.func = [zeros(1, dim.nz*dim.Np), zeros(1, dim.nd*dim.Np), zeros(1, dim.nu*dim.Np), ...
              ones(1, dim.nq*dim.Np), params.lambda*ones(1, dim.np*dim.Np)]';
    
    obj.intcon = dim.nz*dim.Np + 1:(dim.nz + dim.nd)*dim.Np;
    obj.lb = [-inf*ones(dim.nz*dim.Np,1); zeros(length(obj.intcon),1); -inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];
    obj.ub = [inf*ones(dim.nz*dim.Np,1); ones(length(obj.intcon),1); inf*ones((dim.nu + dim.nq + dim.np)*dim.Np,1)];
end