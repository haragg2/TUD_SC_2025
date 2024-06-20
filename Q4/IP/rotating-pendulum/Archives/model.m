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
T2 = 1/2 * I2 * theta2_dot^2 + 1/2 * m2 * c2^2 * theta2_dot^2 + ...
     1/2 * m2 * (l1^2 * theta1_dot^2 + 2 * c2 * l1 * cos(theta2 - theta1) * theta1_dot * theta2_dot);
T = T1 + T2;

% Potential energy terms
U1 = -m1 * g * c1 * cos(theta1);
U2 = -m2 * g * (l1 * cos(theta1) + c2 * cos(theta2));
U = U1 + U2;

% Frictional terms
D1 = 1/2 * b1 * theta1_dot^2;
D2 = 1/2 * b2 * (theta2_dot - theta1_dot)^2;

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

% % Euler-Lagrange equations with friction and external torque
% eq1 = simplify(tau_ext == dL_dtheta1_dot_dt - dL_dtheta1 + dD1_dtheta1_dot + dD2_dtheta1_dot);
% eq2 = simplify(0 == dL_dtheta2_dot_dt - dL_dtheta2 + dD2_dtheta2_dot);
% 
% % Display the Euler-Lagrange equations
% disp('Euler-Lagrange equation for theta1:');
% disp(eq1);
% disp('Euler-Lagrange equation for theta2:');
% disp(eq2);


%% Model
A = b + b1 + b2;
B = I1 + m1*c1^2 + m2*l1^2;
C = m1*c1 + m2*l1;
D = l1*c2*m2;
E = b2;
F = a;
G = m2*c2^2 + I2;
H = m2*c2;


% A = 7.29559;                                                                                                                                                                                                                               
% B = 0.0823764;                                                                                                                                                                                                                               
% C = -0.00336043;                                                                                                                                                                                                                
% D = 0.0142557;                                                                                                                                                                                                                              
% E = 0.00861626;                                                                                                                                                                                                                          
% F = 41.7247;                                                                                                                                                                                                                             
% G = 0.0117936;                                                                                                                                                                                                                                
% H = 0.130666;

M_l = [B, D*cos(theta1(t) - theta2(t));
    D*cos(theta1(t) - theta2(t)), G];

C_l = [A*theta1_dot(t) + D*sin(theta1(t) - theta2(t))*theta2_dot(t)^2 - E*theta2_dot(t); 
       -E*(theta1_dot(t) - theta2_dot(t)) - D*sin(theta1(t) - theta2(t))*theta1_dot(t)^2];

G_l = [C*sin(theta1(t));
        H*sin(theta2(t))];

B_l = [F*u; 0];

theta_ddot = symmatrix(inv(M_l)*(B_l - C_l - G_l*g));

I1 = 0.074;
I2 = 0.000160602;
m1 = 0.125;
m2 = 0.0869285;
c1 = 0.04; % could be negative
c2 = 0.0466771;
l1 = 0.1;
g = 9.81;
b1 = 4.8;
b2 = 6.5799e-05;
a = 50;
b = 8;

eq1 = theta1_ddot == theta_ddot(1, 1);
eq2 = theta2_ddot == theta_ddot(2, 1);

eq1 = subs(eq1);
eq2 = subs(eq2);

[V, S] = odeToVectorField(eq2, eq1);
MF = matlabFunction(V, 'vars', {'t','Y', 'u'});

% u_data = load('input_data_chirp_0.25_1_45.mat');
% theta1_data = load('theta1_chirp.mat');
% theta2_data = load('theta2_chirp.mat');
% all_data = load('data_all.mat');
all_data = load('all_0.25_1_45.mat');

% input_voltage = u_data.input_data.Data;
% input_voltage = all_data.volt.Data(:);
% time = all_data.volt.Time;
% h = time(2) - time(1);
% t_end = time(end);

th_1 = train_datad.y(:, 1);
th_2 = train_datad.y(:, 2);
input_voltage = train_datad.InputData;
h = 0.01;
t_end = 20;
time = 0:h:t_end;

theta0 = [-th_1(1); 0; -th_1(1)-th_2(1); 0];

Theta_eval = [theta0'];

for i=1:1:size(time,2)
    t = time(i);
    u = input_voltage(i);

    tspan = [t, t+h];
    [~, theta] = ode45(@(t, Y)MF(t, Y, u), tspan, theta0);
    theta0 = theta(end, :)';

    Theta_eval = [Theta_eval; theta0'];

end

close all;
figure;
plot(-th_1);
hold on;
plot(Theta_eval(:,1));
hold off;

figure;
plot(-th_1-th_2);
hold on;
plot(Theta_eval(:,3));