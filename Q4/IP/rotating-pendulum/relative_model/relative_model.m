% Define symbolic variables
syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u theta1(t) theta2(t)

% Derivatives of theta1 and theta2
theta1_dot = diff(theta1, t);
theta2_dot = diff(theta2, t);

% Derivatives of theta1_dot and theta2_dot
theta1_ddot = diff(theta1_dot, t);
theta2_ddot = diff(theta2_dot, t);

% Kinetic energy terms
T1 = 1/2 * I1 * theta1_dot^2 + 1/2 * m1 * c1^2 * theta1_dot^2;
T2 = 1/2 * I2 * (theta1_dot + theta2_dot)^2 + ...
     1/2 * m2 * (c2^2 * (theta1_dot + theta2_dot)^2 + l1^2 * theta1_dot^2 + 2 * c2 * l1 * cos(theta2) * theta1_dot * (theta1_dot + theta2_dot));
T = T1 + T2;

% Potential energy terms
U1 = -m1 * g * c1 * cos(theta1);
U2 = -m2 * g * (l1 * cos(theta1) + c2 * cos(theta1 + theta2));
U = U1 + U2;

% Frictional terms
D1 = 1/2 * b1 * theta1_dot^2;
D2 = 1/2 * b2 * theta2_dot^2;

% Total Lagrangian
L = T - U;

% External torque acting on the first link
tau_ext = a*u - b*theta1_dot;

% Euler-Lagrange equations
% dL/d(theta1) and dL/d(theta2)
dL_dtheta1 = diff(L, theta1);
dL_dtheta2 = diff(L, theta2);

% dL/d(theta1_dot) and dL/d(theta2_dot)
dL_dtheta1_dot = diff(L, theta1_dot);
dL_dtheta2_dot = diff(L, theta2_dot);

% Time derivatives of dL/d(theta1_dot) and dL/d(theta2_dot)
dL_dtheta1_dot_dt = diff(dL_dtheta1_dot, t);
dL_dtheta2_dot_dt = diff(dL_dtheta2_dot, t);

% Rayleigh dissipation function terms (friction)
dD1_dtheta1_dot = diff(D1, theta1_dot);
dD2_dtheta2_dot = diff(D2, theta2_dot);
dD2_dtheta1_dot = diff(D2, theta1_dot);

% Euler-Lagrange equations with friction and external torque
eq1 = simplify(tau_ext == dL_dtheta1_dot_dt - dL_dtheta1 + dD1_dtheta1_dot + dD2_dtheta1_dot);
eq2 = simplify(0 == dL_dtheta2_dot_dt - dL_dtheta2 + dD2_dtheta2_dot);

% Display the Euler-Lagrange equations
disp('Euler-Lagrange equation for theta1:');
disp(eq1);
disp('Euler-Lagrange equation for theta2:');
disp(eq2);


%% Model
syms A B C D E F

% A = I1 + m1*c1^2
% C = l1
% D = m1*c1 
% B = I2 + m2*c2^2;
% E = m2*c2;
% F = b + b1;

M_l = [A + m2*C^2 + B + 2*m2*C*c2*cos(theta2), B + m2*C*c2*cos(theta2);
       B + m2*C*c2*cos(theta2), B];

C_l = [F - 2*m2*C*c2*theta2_dot*sin(theta2), -m2*C*c2*theta2_dot*sin(theta2); 
       m2*C*c2*theta1_dot*sin(theta2), b2];

G_l = [(D + m2*C)*sin(theta1) + E*sin(theta1+theta2);
       E*sin(theta1+theta2)];

B_l = [a*u; 0];


theta_ddot = symmatrix(inv(M_l)*(B_l - C_l*[theta1_dot; theta2_dot] - G_l * g));

%% - Simulate non-linear sysytem 

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

eq1 = theta1_ddot == theta_ddot(1, 1);
eq2 = theta2_ddot == theta_ddot(2, 1);

eq1 = subs(eq1);
eq2 = subs(eq2);

[V, S] = odeToVectorField(eq2, eq1);
MF = matlabFunction(V, 'vars', {'t','Y', 'u'});

all_data = load('all_0.25_1_02.mat');
input_voltage = all_data.volt.Data(:);
time = all_data.volt.Time;
h = time(2) - time(1);
t_end = time(end);

th_1 = all_data.th1.Data(:);
th_2 = all_data.th2.Data(:);

theta0 = [th_1(1); 0; th_2(1); 0];

Theta_eval = [theta0'];
for i=1:1:size(time)
    t = time(i);
    u = input_voltage(i);

    tspan = [t, t+h];
    [~, theta] = ode45(@(t, Y)MF(t, Y, u), tspan, theta0);
    theta0 = theta(end, :)';

    Theta_eval = [Theta_eval; theta0'];

end

close all;
figure;
plot(th_1);
hold on;
plot(Theta_eval(:,1));
hold off;

figure;
plot(th_2);
hold on;
plot(Theta_eval(:,3));


%% Linearize the system

vars = [theta1; theta1_dot; theta2; theta2_dot];
theta_ddot = inv(M_l) * (B_l - C_l *[theta1_dot; theta2_dot] - G_l * g);

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
theta1 = 0;
theta2 = -pi;

A = 0.02;
B = I2 + m2*c2^2;
C = 0.0954944;
D = 0.015;
E = m2*c2;
F = 29.9204;

sys_A = eval(subs(sym_A));
sys_B = eval(subs(sym_B));
sys_C = eye(4);

%% Simulate linear system

u_data = data.volt.Data(:);
theta1 = detrend(data.th1.Data(:));
% theta1 = idfilt(data.th1.Data, 2, 0.003, 'high');
theta2 = detrend(-data.th2.Data(:));
% theta2 = idfilt(-data.th2.Data, 2, 0.003, 'high');
t = data.volt.Time(:);
y = lsim(sys_A, sys_B, sys_C, [], u_data, t);

vaf_value_1 = (1 - var(theta1 - y(:,1)) / var(theta1)) * 100;
vaf_value_2 = (1 - var(theta2 - y(:,3)) / var(theta2)) * 100;

% Plotting theta1 and y(:,1) with labels, title, and legend
figure;
hold on;
plot(theta1);
plot(y(:,1));
xlabel('Time');
ylabel('Value');
title('Comparison of theta_{1}', sprintf('VAF: %.2f%%', vaf_value_1));
legend('theta1_{m}', 'theta1_{s}');
hold off;

% Plotting theta2 and y(:,3) with labels, title, and legend
figure;
hold on;
plot(theta2);
plot(y(:,3));
xlabel('Time');
ylabel('Value');
title('Comparison of theta_{2}', sprintf('VAF: %.2f%%', vaf_value_2));
legend('theta2_{m}', 'theta2_{s}');
hold off;