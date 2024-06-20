
%all_data = load('all_0.25_1_60_02.mat');

% u_struct = load('input_data_chirp_0.25_1_45');
% u_data = u_struct.input_data.Data;
u_data = volt.Data(:);

% theta1_struct = load('theta1_chirp');
% theta1_data = theta1_struct.theta1_data.Data;
theta1 = th1.Data;

% theta2_struct = load('theta2_chirp');
% theta2_data = theta2_struct.theta2_data.Data;
% theta2 = pi + theta2_data - theta1_data;
theta2 = th1.Data - th2.Data;
h = volt.Time(2) - volt.Time(1);

z = iddata([theta1, theta2], u_data, h, 'Name', 'Rot_Pendulum');
z.InputName = 'Voltage';
z.InputUnit =  'V';
z.OutputName = {'Theta 1 Angular position', 'Theta 2 Angular position'};
z.OutputUnit = {'rad', 'rad'};
z.Tstart = 0;
z.TimeUnit = 's';

figure('Name', [z.Name ': Voltage input -> Theta 1 Angular position']);
plot(z(:, 1, 1));   % Plot first input-output pair (Voltage -> Angular position).
figure('Name', [z.Name ': Voltage input -> Theta 2 Angular position']);
plot(z(:, 2, 1));   % Plot second input-output pair (Voltage -> Angular velocity).


% A = b + b1 + b2;
% B = I1 + m1*c1^2 + m2*l1^2;
% C = m1*c1 + m2*l1;
% D = l1*c2*m2;
% E = b2;
% F = a;
% G = m2*c2^2 + I2;
% H = m2*c2;

FileName      = 'rot_pendulum_inertia';       % File describing the model structure.
Order         = [2 1 4];           % Model orders [ny nu nx].
Parameters    = [5.57477, 0.0634785, 0.0129319, 0.0104032, 0.00817803, 31.8741, 0.00894318, 0.1014];         % Initial parameters. Np = 6.
InitialStates = [3.17; 0; 2.1096; 0];            % Initial initial states.
Ts            = 0;                                  % Time-continuous system.

nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
                'Name', 'Rot_Pendulum');
set(nlgr, 'InputName', 'Voltage', 'InputUnit', 'V',               ...
          'OutputName', {'Theta 1 Angular position', 'Theta 2 Angular position'}, ...
          'OutputUnit', {'rad', 'rad'},                         ...
          'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'Theta 1 Angular position' 'Theta 1 Angular velocity' ...
                'Theta 2 Angular position' 'Theta 2 Angular velocity'});
nlgr = setinit(nlgr, 'Unit', {'rad' 'rad/s' 'rad' 'rad/s'});
% nlgr = setpar(nlgr, 'Name', {'K_e' 'm1' 'm2' 'l1' 'l2', 'g'});
% nlgr = setpar(nlgr, 'Unit', {'Vs/rad' 'kg' 'kg' 'm' 'm' 'm/s^2'});

nlgr.SimulationOptions.AbsTol = 1e-3;
nlgr.SimulationOptions.RelTol = 1e-3;
compare(z, nlgr);

nlgr = setinit(nlgr, 'Fixed', {false false false false}); % Estimate the initial states.
opt = nlgreyestOptions('Display', 'on');
opt.SearchOptions.MaxIterations = 30;
nlgr = nlgreyest(z, nlgr, opt);

nlgr.Report
fprintf('\n\nThe search termination condition:\n')
nlgr.Report.Termination

% compare(z, nlgr);