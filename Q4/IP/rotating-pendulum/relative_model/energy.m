% Define symbolic variables
syms I1 I2 m1 m2 c1 c2 l1 g b1 b2 a b u theta1(t) theta2(t) A B C D

% Derivatives of theta1 and theta2
theta1_dot = diff(theta1, t);
theta2_dot = diff(theta2, t);

% Derivatives of theta1_dot and theta2_dot
theta1_ddot = diff(theta1_dot, t);
theta2_ddot = diff(theta2_dot, t);

% Kinetic energy terms
T1 = 1/2 * theta1_dot^2*A;
T2 = 1/2 * I2 * (theta1_dot + theta2_dot)^2 + ...
     1/2 * m2 * (c2^2 * (theta1_dot + theta2_dot)^2 + l1^2 * theta1_dot^2 + 2 * c2 * l1 * cos(theta2) * theta1_dot * (theta1_dot + theta2_dot));
T = T1 + T2;

% Potential energy terms
U1 = -D * g * cos(theta1);
U2 = -m2 * g * (l1 * cos(theta1) + c2 * cos(theta1 + theta2));
U = U1 + U2;

a = 187.247;
g = 9.81;
I2 = 0.000110802;
m2 = 0.05914;
c2 = 0.0563237;
b2 = 5.65455e-05;

A = 0.02;
B = I2 + m2*c2^2;
l1 = 0.0954944;
D = 0.015;
E = m2*c2;
F = 29.9204;

theta1 = pi;
theta2 = 0;
theta1_dot = 0;

theta2_dot = 0;


T = subs(T)
U = eval(subs(U))