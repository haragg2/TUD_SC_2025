syms theta(t) theta_p(t) m_w m_p g l R beta J Tau

theta_d = diff(theta);
theta_pd = diff(theta_p);

theta_dd = diff(theta_d);
theta_pdd = diff(theta_pd);

eq1 = theta_dd(t) == ((m_p * l^2 + m_p * R * l * cos(theta_p(t)) * cos(beta)) ...
    * Tau + (m_p * l^2) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2 - (m_p * R * l * cos(theta_p(t)) ...
    * cos(beta)) * m_p * g * l * sin(theta_p(t)) - (m_p * l^2) * ((m_w + m_p) * g * R * ...
    sin(beta))) / ((J + ((m_w + m_p) * R^2)) * (m_p * l^2) - (m_p * R * l * cos(theta_p(t)) * cos(beta))^2);


eq2 = theta_pdd(t) == ((-(J + ((m_w + m_p) * R^2)) - (m_p * R * l * cos(theta_p(t)) * cos(beta))) ...
    * Tau + (J + ((m_w + m_p) * R^2)) * m_p * g * l * sin(theta_p(t)) - (m_p * R * l * cos(theta_p(t)) ...
    * cos(beta)) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2 + (m_p * R * l * cos(theta_p(t)) ...
    * cos(beta)) * ((m_w + m_p) * g * R * sin(beta)))  / ((J + ((m_w + m_p) * R^2)) ...
    * (m_p * l^2) - (m_p * R * l * cos(theta_p(t)) * cos(beta))^2);


% Model Parameters
m_w = 1.5;          % kg; mass of wheel
m_p = 6.16;         % kg; mass of pendulum
l = 0.46;           % m; distance between COG of pendulum and center of wheel
R = 0.054;          % m; radius of wheel
g = 9.81;           % m/s^2; gravity
J = 0.1;            % kgm^2; wheel centroidal inertia
beta = 10 * (pi/180);
Tau = 0;
% Tau = (m_w + m_p) * g * R * sin(beta);

eq1 = subs(eq1);
eq2 = subs(eq2);

[V,S] = odeToVectorField(eq2, eq1);
M = matlabFunction(V,'vars',{'t','Y'});

theta_p_eq = asin((m_w + m_p) * R * sin(beta) / (m_p * l));

t_end = 10;

initCond = [0.1 0 theta_p_eq 0];
sols = ode45(M, [0 t_end], initCond);

figure;
plot(sols.x, sols.y)
legend('\theta_p','d\theta_p/dt','\theta','d\theta/dt')
title('Solutions of State Variables')
xlabel('Time (s)')
ylabel('Solutions (rad or rad/s)')

xw = @(t) (R * cos(beta) * deval(sols, t, 1));
yw = @(t) (R * sin(beta) * deval(sols, t, 1));

xp = @(t) ((R * cos(beta) * deval(sols, t, 1)) + (l * sin(deval(sols, t, 3))));
yp = @(t) ((R * sin(beta) * deval(sols, t, 1)) + (l * cos(deval(sols, t, 3))));

figure;
fanimator(@(t) plot(xp(t), yp(t),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 t_end]);
axis equal;
hold on;
fanimator(@(t) plot([xw(1) xw(t)],[yw(1) yw(t)],'b-'), 'AnimationRange', [0 t_end]);
% fanimator(@(t) plot([xw(t) xw(t) - l/2],[yw(t) yw(t) + l/2],'k-'), 'AnimationRange', [0 t_end]);
fanimator(@(t) plot([xw(t) xp(t)],[yw(t) yp(t)],'k-'), 'AnimationRange', [0 t_end]);
fanimator(@(t) plot(xw(t), yw(t),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 t_end]);
fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)), 'AnimationRange', [0 t_end]);
hold off;

