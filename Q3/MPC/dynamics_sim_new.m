syms theta(t) theta_p(t) m_w m_p g l R beta J Tau

theta_d = diff(theta);
theta_pd = diff(theta_p);

theta_dd = diff(theta_d);
theta_pdd = diff(theta_pd);

eq1 = theta_dd(t) == (((m_p * l^2) + (m_p * R * l * cos(theta_p(t) + beta))) * Tau ...
    + ((m_p * l^2) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2) ...
    - (m_p * R * l * cos(theta_p(t) + beta) * m_p * g * l * sin(theta_p(t))) ...
    - (m_p * l^2) * ((m_w + m_p) * g * R * sin(beta))) / (((J ...
    + ((m_w + m_p) * R^2)) * (m_p * l^2)) - ((m_p * R * l * cos(theta_p(t) + beta))^2));


eq2 = theta_pdd(t) == (((-(J + ((m_w + m_p) * R^2)) - (m_p * R * l * cos(theta_p(t) + beta))) * Tau) ...
    + ((J + ((m_w + m_p) * R^2)) * m_p * g * l * sin(theta_p(t))) ...
    - ((m_p * R * l * cos(theta_p(t) + beta)) * m_p * R * l * sin(theta_p(t)) * theta_pd(t)^2) ...
    + ((m_p * R * l * cos(theta_p(t) + beta)) * ((m_w + m_p) * g * R * sin(beta))))  / (((J ...
    + ((m_w + m_p) * R^2)) * (m_p * l^2)) - ((m_p * R * l * cos(theta_p(t) + beta))^2));


% Model Parameters
m_w = 1.5;          % kg; mass of wheel
m_p = 6.16;         % kg; mass of pendulum
l = 0.46;           % m; distance between COG of pendulum and center of wheel
R = 0.054;          % m; radius of wheel
g = 9.81;           % m/s^2; gravity
J = 0.1;            % kgm^2; wheel centroidal inertia
beta = 10 * (pi/180);
% Tau = 0;
% Tau = (m_w + m_p) * g * R * sin(beta);
% 
% eq1 = subs(eq1);
% eq2 = subs(eq2);
% 
% [V,S] = odeToVectorField(eq2, eq1);
% M = matlabFunction(V,'vars',{'t','Y', 'Tau'});
% 
theta_p_eq = asin((m_w + m_p) * R * sin(beta) / (m_p * l));

% t_end = 10;

% initCond = [0 0 theta_p_eq 0];
% initCond = [16.9, 3.89, theta_p_eq-0.13, -0.09];

x00 = [16.9, 3.89, -0.11, -0.09]';
x00 = [0, 0, 0, 0]';
xx = zeros(length(sys_dis.A(:,1)), size(T, 2));
xx(:,1) = x00;

for k = 1:1:size(T, 2)-1
    x00 = xx(:, k);
    tspan = [T(k), T(k+1)];
    Tau_eq = (m_w + m_p) * g * R * sin(beta);
    xx(:,k+1) = NLSystemDynamics(x00, tspan, 0);
end
xx(3,:) = xx(3,:) + theta_p_eq;

xw = @(tt) (r * cos(beta) * xx(1, (uint16(tt/0.05)+1)));
yw = @(tt) (r * sin(beta) * xx(1, (uint16(tt/0.05)+1)));

xp = @(tt) ((r * cos(beta) * xx(1, (uint16(tt/0.05)+1)) + (l * sin(xx(3, (uint16(tt/0.05)+1))))));
yp = @(tt) ((r * sin(beta) * xx(1, (uint16(tt/0.05)+1)) + (l * cos(xx(3, (uint16(tt/0.05)+1))))));


% T=0:0.05:10;
% ts = [0];
% xs = [initCond];

% for i=1:1:size(T, 2)-1
%     tspan = [T(i), T(i+1)];
% 
%     Tau = (m_w + m_p) * g * R * sin(beta);
% 
%     [t_step, x_step] = ode45(@(t, Y)M(t, Y, Tau), tspan, initCond, options);
%     initCond = x_step(end, :);
% 
%     ts = [ts; t_step(2:end)]; % append new time points, avoid duplication
%     xs = [xs; x_step(2:end, :)]; % append new results, avoid duplication
% end

% xs = xs(mod(ts, 0.05) == 0, :);

% xw = @(tt) (r * cos(beta) * xs((uint16(tt/0.05)+1), 1));
% yw = @(tt) (r * sin(beta) * xs((uint16(tt/0.05)+1), 1));
% 
% xp = @(tt) ((r * cos(beta) * xs((uint16(tt/0.05)+1), 1) + (l * sin(xs((uint16(tt/0.05)+1), 3)))));
% yp = @(tt) ((r * sin(beta) * xs((uint16(tt/0.05)+1), 1) + (l * cos(xs((uint16(tt/0.05)+1), 3)))));

figure;
axis equal;
hold on;
fanimator(@(tt) plot(xp(tt), yp(tt),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(Ts) xw(tt)],[yw(Ts) yw(tt)],'b-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot([xw(tt) xp(tt)],[yw(tt) yp(tt)],'k-'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) plot(xw(tt), yw(tt),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 T_sim],'FrameRate',20);
fanimator(@(tt) text(-0.3,0.3,"Timer: "+ num2str(tt, 3)), 'AnimationRange', [0 T_sim],'FrameRate',20);
hold off;

% figure;
% plot(sols.x, sols.y)
% legend('\theta_p','d\theta_p/dt','\theta','d\theta/dt')
% title('Solutions of State Variables')
% xlabel('Time (s)')
% ylabel('Solutions (rad or rad/s)')
% 
% xw = @(t) (R * cos(beta) * deval(sols, t, 1));
% yw = @(t) (R * sin(beta) * deval(sols, t, 1));
% 
% xp = @(t) ((R * cos(beta) * deval(sols, t, 1)) + (l * sin(deval(sols, t, 3))));
% yp = @(t) ((R * sin(beta) * deval(sols, t, 1)) + (l * cos(deval(sols, t, 3))));
% 
% figure;
% fanimator(@(t) plot(xp(t), yp(t),'ro','MarkerSize', 10,'MarkerFaceColor','r'), 'AnimationRange', [0 t_end]);
% axis equal;
% hold on;
% fanimator(@(t) plot([xw(1) xw(t)],[yw(1) yw(t)],'b-'), 'AnimationRange', [0 t_end]);
% % fanimator(@(t) plot([xw(t) xw(t) - l/2],[yw(t) yw(t) + l/2],'k-'), 'AnimationRange', [0 t_end]);
% fanimator(@(t) plot([xw(t) xp(t)],[yw(t) yp(t)],'k-'), 'AnimationRange', [0 t_end]);
% fanimator(@(t) plot(xw(t), yw(t),'go','MarkerSize', 10,'MarkerFaceColor','g'), 'AnimationRange', [0 t_end]);
% fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)), 'AnimationRange', [0 t_end]);
% hold off;
% 

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