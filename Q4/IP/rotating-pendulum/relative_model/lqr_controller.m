
%% Linearize the system
clc;
% Linearization point
theta1 = 0;
theta2 = 0;

%[A, B, C, D] = LinearizeSystem(theta1, theta2);
C = [1 0 0 0; 0 0 1 0];

load("LinearSys_pi_0.mat");
A = sys_A;
B = sys_B;
D = [];


sys_ct = ss(A, B, C, D);

%h = 0.1;
%% Discretize the system
clc;
sys_dis = c2d(sys_ct, h);

%pi_0
% R = 0.1;
% 
% Q = [40 0 0 0;
%      0 10 0 0;
%      0 0 70 0;
%      0 0 0 10];

%0_-pi
R = 0.1;

Q = [40 0 0 0;
     0 10 0 0;
     0 0 70 0;
     0 0 0 10];

% LQI
for i=0.01:0.01:1.8
    Q_lqi = blkdiag(Q, i*eye(2));
    try
        [K_lqr,~,~] = lqi(sys_dis,Q_lqi,R);
    catch
        i;
    end
end


% K_lqr = dlqr(sys_dis.A, sys_dis.B, Q, R);
% K_lqr = [K_lqr, 0, 0];
sys_cl = ss(sys_dis.A-sys_dis.B*K_lqr(1:4), sys_dis.B, sys_dis.C, sys_dis.D, h);

% Get the model function
%model = getModel();

% tspan = 0:h:10;
% u = zeros(length(tspan),1);
% x0 = [.2;0;.1;0];
% lsim(sys_cl, u', tspan, x0)
%% Simulation
% t_end = 5;
% T = 0:h:t_end;
% 
% usim = zeros(1, size(T, 2)-1);
% t = zeros(size(T, 2)-1, 1);
% 
% x0 = [0, 0, 0.32, 0]';  % initial condition
% x = zeros(4, size(T, 2));
% x(:,1) = x0;
% 
% x_eq = [pi; 0; 0; 0];
% 
% for k = 1:1:size(T, 2)-1
%     t(k) = (k-1) * h;
%     if ( k > 1 && (floor(t(k)) - floor(t(k-1))) == 1 )
%         fprintf('t = %d sec \n', floor(t(k)));
%     end
% 
%     % Get the control input from the LQR controller
%     usim(:,k) = -K_lqr * x(:,k);
% 
%     tspan = [T(k), T(k+1)];
%     % x(:,k+1) = NLSystemDynamics(x_eq, x(:,k), tspan, usim(:,k));
%     x(:,k+1) = NLSystemDynamics(model, x_eq, x(:,k), tspan, usim(:,k));
% 
% end
% 
% 
% figure;
% sgtitle("Regulation LQR");
% subplot(5, 1, 1);
% hold on;
% grid on;
% hplots(1) = stairs(x(1,1:end-1), 'LineWidth', 1.5);
% % hplots(2) = stairs(x_LQ(1:end-1,1), 'LineWidth', 1.5);
% ylabel('$\theta_{1}$');
% hold off;
% subplot(5, 1, 2);
% hold on;
% grid on;
% hplots(1) = stairs(x(2,1:end-1), 'LineWidth', 1.5);
% % hplots(2) = stairs(x_LQ(1:end-1,2), 'LineWidth', 1.5);
% ylabel('$$\dot{\theta_{1}}$$');
% hold off;
% subplot(5, 1, 3);
% hold on;
% grid on;
% hplots(1) = stairs(x(3,1:end-1), 'LineWidth', 1.5);
% % hplots(2) = stairs(x_LQ(1:end-1,3), 'LineWidth', 1.5);
% ylabel('$\theta_{2}$');
% hold off;
% subplot(5, 1, 4);
% hold on;
% grid on;
% hplots(1) = stairs(x(4,1:end-1), 'LineWidth', 1.5);
% % hplots(2) = stairs(x_LQ(1:end-1,4), 'LineWidth', 1.5);
% ylabel('$$\dot{\theta_{2}}$$');
% hold off;
% subplot(5, 1, 5);
% hold on;
% grid on;
% hplots(1) = stairs(usim, 'LineWidth', 1.5);
% % hplots(1) = stairs(u_LQ(1:end-1), 'LineWidth', 1.5);
% hold off;
% xlabel('Time Step $k$');
% ylabel('$u$');

%% Functions


% function [X] = NLSystemDynamics(x_eq, x0, tspan, input_voltage)
% 
%     syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u A B C D E F theta1(t) theta2(t)
% 
%     % A = I1 + m1*c1^2
%     % C = l1
%     % D = m1*c1 
%     % B = I2 + m2*c2^2;
%     % E = m2*c2;
%     % F = b + b1;
% 
%     theta1_dot = diff(theta1, t);
%     theta2_dot = diff(theta2, t);
% 
%     % Derivatives of theta1_dot and theta2_dot
%     theta1_ddot = diff(theta1_dot, t);
%     theta2_ddot = diff(theta2_dot, t);
% 
%     M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
%            B + m2*C*c2*cos(theta2), B];
% 
%     C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
%            m2*C*c2*theta1_dot*sin(theta2), b2];
% 
%     G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
%            E*sin(theta1+theta2)];
% 
%     B_l = [a*u; 0];
% 
% 
%     theta_ddot = symmatrix(inv(M_l)*(B_l - C_l*[theta1_dot; theta2_dot] - G_l * g));
% 
%     % Identified parameters - copy paste
%     a = 187.247;
%     g = 9.81;
%     I2 = 0.000110802;
%     m2 = 0.05914;
%     c2 = 0.0563237;
%     b2 = 5.65455e-05;
% 
%     A = 0.02;
%     B = I2 + m2*c2^2;
%     C = 0.0954944;
%     D = 0.015;
%     E = m2*c2;
%     F = 29.9204;
% 
%     eq1 = theta1_ddot == theta_ddot(1, 1);
%     eq2 = theta2_ddot == theta_ddot(2, 1);
% 
%     eq1 = subs(eq1);
%     eq2 = subs(eq2);
% 
%     [V, S] = odeToVectorField(eq2, eq1);
%     MF = matlabFunction(V, 'vars', {'t','Y', 'u'});
% 
%     x0 = x0';
%     x0 = x0 + x_eq';
%     u = input_voltage;
%     [~, x_step] = ode45(@(t, Y)MF(t, Y, u), tspan, x0);
% 
%     X = x_step(end, :)';
%     X = X - x_eq;
% end


function [sys_A, sys_B, sys_C, sys_D] = LinearizeSystem(th1, th2)

    syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u A B C D E F theta1(t) theta2(t)
    
    % A = I1 + m1*c1^2
    % C = l1
    % D = m1*c1 
    % B = I2 + m2*c2^2;
    % E = m2*c2;
    % F = b + b1;


    theta1_dot = diff(theta1, t);
    theta2_dot = diff(theta2, t);

    M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
           B + m2*C*c2*cos(theta2), B];
    
    C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
           m2*C*c2*theta1_dot*sin(theta2), b2];
    
    G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
           E*sin(theta1+theta2)];
    
    B_l = [a*u; 0];
    
    theta_ddot = inv(M_l) * (B_l - C_l *[theta1_dot; theta2_dot] - G_l * g);
    vars = [theta1; theta1_dot; theta2; theta2_dot];
    
    % Compute the Jacobian
    A_jac = symmatrix(jacobian(theta_ddot, vars));
    
    sym_A = [jacobian(theta1_dot,vars); A_jac(1,:);
             jacobian(theta2_dot,vars); A_jac(2,:);];
    %sym_A = symmatrix2sym(sym_A);
    
    B_jac = symmatrix(jacobian(theta_ddot, u));
    
    sym_B = [jacobian(theta1_dot,u); B_jac(1,:);
             jacobian(theta2_dot,u); B_jac(2,:);];
    %sym_B = symmatrix2sym(sym_B);
    
    theta1 = th1;
    theta2 = th2;

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

    % Evaluate the matrices
    sys_A = eval(subs(sym_A));
    sys_B = eval(subs(sym_B));
    sys_C = [1, 0, 0, 0;
             0, 0, 1, 0];
    sys_D = [];
end


function model = getModel()

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
    
    theta_ddot = symmatrix(inv(M_l)*(B_l - C_l*[theta1_dot; theta2_dot] - G_l * g));
    
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
    
    eq1 = theta1_ddot == theta_ddot(1, 1);
    eq2 = theta2_ddot == theta_ddot(2, 1);
    
    eq1 = subs(eq1);
    eq2 = subs(eq2);
    
    [V, S] = odeToVectorField(eq2, eq1);
    MF = matlabFunction(V, 'vars', {'t','Y', 'u'});

    model = MF;
end



function [X] = NLSystemDynamics(model, x_eq, x0, tspan, input_voltage)
    
    x0 = x0';
    x0 = x0 + x_eq';
    u = input_voltage;
    [~, x_step] = ode45(@(t, Y) model(t, Y, u), tspan, x0);
    
    X = x_step(end, :)'
    X = X - x_eq;
end
