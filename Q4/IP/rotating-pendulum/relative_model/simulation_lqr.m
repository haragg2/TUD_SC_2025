close all;
clear all;
clc;

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Linearize the system

% Linearization point
theta1_eq = pi;
theta2_eq = 0;

[A, B, C, D] = LinearizeSystem(theta1_eq, theta2_eq);
sys_ct = ss(A, B, C, D);

h = 0.01;

%% Discretize the system
sys_dis = c2d(sys_ct, h);

Q = [40 0 0 0;
     0 10 0 0;
     0 0 70 0;
     0 0 0 10];

R = 0.1;

K_lqr = dlqr(sys_dis.A, sys_dis.B, Q, R);
sys_cl = ss(sys_dis.A - sys_dis.B*K_lqr, sys_dis.B, sys_dis.C, sys_dis.D, h);

% Get the model function
model = getModel();

% all_data = load('../data/all_0.25_1_02.mat');
% usim = 5*all_data.volt.Data(:);
% T = all_data.volt.Time;
% theta1_0 = detrend(all_data.th1.Data(:));
% theta2_0 = detrend(all_data.th2.Data(:));
% t_end = T(end);

% x0 = [0.1, 0, 0, 0]';  % initial condition
% x = zeros(4, size(T, 2));
% x(:,1) = x0;

t_end = 5;
T = 0:h:t_end;

usim = zeros(1, size(T, 2)-1);
t = zeros(size(T, 2)-1, 1);

x0 = [0.4, 0, -0.1, 0]';  % initial condition
x = zeros(4, size(T, 2));
x(:,1) = x0;

x_eq = [theta1_eq; 0; theta2_eq; 0];

for k = 1:1:size(T, 2)-1
    t(k) = (k-1) * h;
    if ( k > 1 && (floor(t(k)) - floor(t(k-1))) == 1 )
        fprintf('t = %d sec \n', floor(t(k)));
    end

    % Get the control input from the LQR controller
    usim(:,k) = -K_lqr * x(:,k);
    usim(:,k) = max(min(usim(:,k), 1), -1);

    tspan = [T(k), T(k+1)];
    % x(:,k+1) = NLSystemDynamics(x_eq, x(:,k), tspan, usim(:,k));
    x(:,k+1) = NLSystemDynamics(model, x_eq, x(:,k), tspan, usim(:, k));

end


figure;
sgtitle("Regulation LQR");
subplot(5, 1, 1);
hold on;
grid on;
hplots(1) = stairs(x(1,1:end-1), 'LineWidth', 1.5);
% hplots(2) = stairs(x_LQ(1:end-1,1), 'LineWidth', 1.5);
ylabel('$\theta_{1}$');
hold off;
subplot(5, 1, 2);
hold on;
grid on;
hplots(1) = stairs(x(2,1:end-1), 'LineWidth', 1.5);
% hplots(2) = stairs(x_LQ(1:end-1,2), 'LineWidth', 1.5);
ylabel('$$\dot{\theta_{1}}$$');
hold off;
subplot(5, 1, 3);
hold on;
grid on;
hplots(1) = stairs(x(3,1:end-1), 'LineWidth', 1.5);
% hplots(2) = stairs(x_LQ(1:end-1,3), 'LineWidth', 1.5);
ylabel('$\theta_{2}$');
hold off;
subplot(5, 1, 4);
hold on;
grid on;
hplots(1) = stairs(x(4,1:end-1), 'LineWidth', 1.5);
% hplots(2) = stairs(x_LQ(1:end-1,4), 'LineWidth', 1.5);
ylabel('$$\dot{\theta_{2}}$$');
hold off;
subplot(5, 1, 5);
hold on;
grid on;
hplots(1) = stairs(usim, 'LineWidth', 1.5);
% hplots(1) = stairs(u_LQ(1:end-1), 'LineWidth', 1.5);
hold off;
xlabel('Time Step $k$');
ylabel('$u$');

%%

setupAnimation(1, x, x_eq, h, t_end);
%% Functions

function [sys_A, sys_B, sys_C, sys_D] = LinearizeSystem(th1, th2)

    % syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u A B C D E F theta1(t) theta2(t)
    % 
    % % A = I1 + m1*c1^2
    % % C = l1
    % % D = m1*c1 
    % % B = I2 + m2*c2^2;
    % % E = m2*c2;
    % % F = b + b1;
    % 
    % theta1_dot = diff(theta1, t);
    % theta2_dot = diff(theta2, t);
    % 
    % M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
    %        B + m2*C*c2*cos(theta2), B];
    % 
    % C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
    %        m2*C*c2*theta1_dot*sin(theta2), b2];
    % 
    % G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
    %        E*sin(theta1+theta2)];
    % 
    % B_l = [a*u; 0];
    % 
    % theta_ddot = inv(M_l) * (B_l - C_l *[theta1_dot; theta2_dot] - G_l * g);

    % syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u A B C D E F theta1(t) theta2(t)
    syms u theta1(t) theta2(t)
    theta1_dot = diff(theta1, t);
    theta2_dot = diff(theta2, t);
    [~, theta_ddot] = getModel();
    vars = [theta1; theta1_dot; theta2; theta2_dot];
    
    % Compute the Jacobian
    A_jac = symmatrix(jacobian(theta_ddot, vars));
    
    sym_A = [jacobian(theta1_dot,vars); A_jac(1,:);
             jacobian(theta2_dot,vars); A_jac(2,:);];
    
    B_jac = symmatrix(jacobian(theta_ddot, u));
    
    sym_B = [jacobian(theta1_dot,u); B_jac(1,:);
             jacobian(theta2_dot,u); B_jac(2,:);];
    
    % Identified parameters
    a = 187.247;
    g = 9.81;
    I2 = 0.000110802;
    m2 = 0.05914;
    c2 = 0.0563237;
    b2 = 5.65455e-05;
    
    A = 0.02;
    B = I2 + m2*c2^2;
    C = 0.0954944;
    D = 0.015;
    E = m2*c2;
    F = 29.9204;

    theta1 = th1;
    theta2 = th2;
    
    % Evaluate the matrices
    sys_A = eval(subs(sym_A));
    sys_B = eval(subs(sym_B));
    sys_C = [1, 0, 0, 0;
             0, 0, 1, 0];
    sys_D = [];
end


function [model, theta_ddot] = getModel()

    syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u A B C D E F theta1(t) theta2(t)
    
    % A = I1 + m1*c1^2
    % C = l1
    % D = m1*c1 
    % B = I2 + m2*c2^2;
    % E = m2*c2;
    % F = b + b1;

    theta1_dot = diff(theta1, t);
    theta2_dot = diff(theta2, t);

    % Derivatives of theta1_dot and theta2_dot
    theta1_ddot = diff(theta1_dot, t);
    theta2_ddot = diff(theta2_dot, t);
    
    M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
           B + m2*C*c2*cos(theta2), B];
    
    C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
           m2*C*c2*theta1_dot*sin(theta2), b2];
    
    G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
           E*sin(theta1+theta2)];
    
    B_l = [a*u; 0];
    
    theta_ddot = inv(M_l)*(B_l - C_l*[theta1_dot; theta2_dot] - G_l * g);
    theta_ddot_s = symmatrix(theta_ddot);
    
    % Identified parameters - copy paste
    a = 187.247;
    g = 9.81;
    I2 = 0.000110802;
    m2 = 0.05914;
    c2 = 0.0563237;
    b2 = 5.65455e-05;
    
    A = 0.02;
    B = I2 + m2*c2^2;
    C = 0.0954944;
    D = 0.015;
    E = m2*c2;
    F = 29.9204;
    
    eq1 = theta1_ddot == theta_ddot_s(1, 1);
    eq2 = theta2_ddot == theta_ddot_s(2, 1);
    
    eq1 = subs(eq1);
    eq2 = subs(eq2);
    
    [V, S] = odeToVectorField(eq2, eq1);
    MF = matlabFunction(V, 'vars', {'t','Y', 'u'});

    model = MF;
end

function [X] = NLSystemDynamics(model, x_eq, x0, tspan, input_voltage)
    
    x0 = x0' + x_eq';
    u = input_voltage;
    [~, x_step] = ode45(@(t, Y) model(t, Y, u), tspan, x0);
    X = x_step(end, :)' - x_eq;
end

function setupAnimation(sm_nl_sys, x, theta_eq, Ts, T_sim)
    l1 = 0.1;
    l2 = 0.1;
    
    % Calculate FPS
    Fps = floor(1/Ts);

    % Add the eqbm pendulum angle for the non-linear system (linear to non-linear translation)
    if sm_nl_sys == 1
        x = x + theta_eq;
    end
    
    % Non-linear equations
    x1 = @(tt) (-l1 * sin(x(1, (floor(tt/Ts)+1))));
    y1 = @(tt) (-l1 * cos(x(1, (floor(tt/Ts)+1))));
    
    x2 = @(tt) (-l2 * sin(x(1, (floor(tt/Ts)+1)) + x(3, (floor(tt/Ts)+1))));
    y2 = @(tt) (-l2 * cos(x(1, (floor(tt/Ts)+1)) + x(3, (floor(tt/Ts)+1))));

    figure;
    xlabel("Horizontal Distance (m)");
    ylabel("Vertical Distance (m)");
    title("Rotated Pendulum");
    xlim([-0.3 0.3]);
    ylim([-0.3 0.3]);

    hold on;
    fanimator(@(tt) plot([0 x1(tt)], [0 y1(tt)],'b-','LineWidth',2), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    fanimator(@(tt) plot([x1(tt) (x1(tt)+x2(tt))], [y1(tt) (y1(tt)+y2(tt))],'b-','LineWidth',2), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    fanimator(@(tt) plot(x1(tt), y1(tt),'ro','MarkerSize', 18,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    fanimator(@(tt) plot((x1(tt)+x2(tt)), (y1(tt)+y2(tt)),'go','MarkerSize', 12,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    fanimator(@(tt) plot(0, 0,'ko','MarkerSize', 5,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    
    % Add the timer
    fanimator(@(tt) text(0,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',Fps);
    hold off;
end