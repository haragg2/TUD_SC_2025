% Define initial parameters
kt = 5;   % Damping coefficient
ke = 0.01;   % Damping coefficient
m1 = 0.2; % Mass of first pendulum
m2 = 0.04; % Mass of second pendulum
l1 = 0.5; % Length of first pendulum
l2 = 0.1; % Length of second pendulum
g = 9.81; % Acceleration due to gravity

% Define initial state [theta1, theta1_dot, theta2, theta2_dot]
x0 = [0.7090; 0; -0.7085; 0]; % Example initial conditions

% Define time span for the simulation
tspan = [0 60]; % Simulate from 0 to 10 seconds

% Define the control input u (e.g., torque applied)
u = input_data.Data; % External voltage

% Function handle to the pendulum dynamics
pendulumDynamics = @(t, x) rot_pendulum(t, x, u, kt, ke, m1, m2, l1, l2, g);

% Solve the system using ode45
[t, x] = ode45(pendulumDynamics, tspan, x0);

% Plot the results
figure;
subplot(2,1,1);
plot(t, x(:,1), 'r-', t, x(:,3), 'b-');
title('Pendulum Angles');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Theta1', 'Theta2');

subplot(2,1,2);
plot(t, x(:,2), 'r-', t, x(:,4), 'b-');
title('Angular Velocities');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('Theta1 dot', 'Theta2 dot');


function dx = rot_pendulum(t, x, u, kt, ke, m1, m2, l1, l2, g, varargin)
    % Output equations.
    y = [x(1);                         
       x(3)                          
      ];

    % State equations.
    dx = [x(2);                        
           (kt*(u(1) - ke*x(2)) - m2*l2*y(2)*cos(x(1) - x(3)) - ...
                m2*l2*x(4)^2*sin(x(1) - x(3)) - ...
                g*(m1 + m2)*sin(x(1))) / ((m1 + m2)*l1); 
           x(4);
           (- m2*l1*y(1)*cos(x(1) - x(3)) - ...
                m2*l1*x(2)^2*sin(x(1) - x(3)) - ...
                g*m2*sin(x(3))) / (m2*l2);    
          ];
end