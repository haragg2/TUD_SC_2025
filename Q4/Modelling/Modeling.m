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

tsim = 0:h:10;
v_init = 20;
X0 = [0; v_init];
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


v_init = 20;
x_init = 0.1;
X0 = [x_init; v_init];

tspan = [tsim(1) 200];
[t_model, X_model] = ode45(@(t, Y) pwa_model(t, Y, params, umax), tspan, X0);
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

%%

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

% % Function to solve the ODE and capture theta
% function dYdt = odefun(t, Y)
%     u = sin(t);  % Example control input, replace with appropriate function
%     [dx, theta] = pwa_model(t, Y, params, u);
%     theta_values = [theta_values; theta];  % Store theta at each step
%     dYdt = dx;
% end


% % Define symbolic variables
% syms x(t) v(t) u(t) r
% 
% theta = 0;
% 
% % Define the differential equations
% ode1 = diff(x, t) == v;
% ode2 = diff(v, t) == ((b/m) * u(t)) / (1 + gamma * r) - (g * sin(theta) * x(t)) - (c/m) * v^2;
% 
% % Convert the system of ODEs to MATLAB function handle
% odes = [ode1; ode2];
% [V, S] = odeToVectorField(odes);
% 
% 
% MF = matlabFunction(V, 'vars', {'t', 'Y', 'u', 'r'});
% 
% h = 0.1;
% tsim = 0:h:10;
% usim = sin(tsim);
% % r_sim = ones(size(tsim));  % Example r(t) function
% v_init = 0;
% r = 1;
% X0 = [0; v_init];
% 
% % Loop through the simulation time
% for i = 1:length(tsim)-1
%     tspan = [tsim(i), tsim(i+1)];
%     u = usim(i);
% 
%     % Define options for ODE solver with events
%     options = odeset('Events', @(t, x) event_function(t, x, vg));
%     % Solve the ODEs
%     [t, X, te, xe, ie] = ode45(@(t, Y) MF(t, Y, u, r), tspan, X0, options);
% 
%     if ~isempty(te)
%         switch ie(end)  % Handle the last event
%             case 1
%                 disp('Event: v > vg');
%             case 2
%                 disp('Event: v <= vg');
%         end
%     end
% 
%     % Update the initial condition for the next step
%     X0 = X(end, :);
% end
% 
% function [value, isterminal, direction] = event_function(t, x, vg)
%     % Define event conditions
%     value = [x(2) - vg, x(2) - vg];  % Detect when v > vg and v <= vg
%     isterminal = [1, 1];  % Stop the integration
%     direction = [1, -1];  % Detect both increasing and decreasing
% end