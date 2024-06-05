
%all_data = load('all_0.25_1_60_02.mat');

% u_struct = load('input_data_chirp_0.25_1_45');
% u_data = u_struct.input_data.Data;
u_temp = volt.Data(:);
u_data = u_temp(215:1442);

% theta1_struct = load('theta1_chirp');
% theta1_data = theta1_struct.theta1_data.Data;
theta1 = th1.Data;

%theta2 = detrend(th2.Data(310:1981)); for th2_0 = 0.5;
theta2 = detrend(th2.Data(215:1442));
h = volt.Time(2) - volt.Time(1);

z = iddata(theta2, u_data, h, 'Name', 'Rot_Pendulum');
z.InputName = 'Voltage';
z.InputUnit =  'V';
% z.OutputName = {'Theta 1 Angular position', 'Theta 2 Angular position'};
% z.OutputUnit = {'rad', 'rad'};
z.Tstart = 0;
z.TimeUnit = 's';

figure('Name', [z.Name ': Voltage input -> Theta 1 Angular position']);
plot(z(:, 1, 1));   % Plot first input-output pair (Voltage -> Angular position).
% figure('Name', [z.Name ': Voltage input -> Theta 2 Angular position']);
% plot(z(:, 2, 1));   % Plot second input-output pair (Voltage -> Angular velocity).


% A = b + b1 + b2;
% B = I1 + m1*c1^2 + m2*l1^2;
% C = m1*c1 + m2*l1;
% D = l1*c2*m2;
% E = b2;
% F = a;
% G = m2*c2^2 + I2;
% H = m2*c2;

FileName      = 'rot_pendulum_theta2';       % File describing the model structure.
Order         = [1 1 2];           % Model orders [ny nu nx].

%Parameters    = [0.054797, 0.0572183, 6.62561e-05,  0.0001]; % t2_0 = 0.3 85% 
Parameters    = [0.05914, 0.0563237, 5.65455e-05,  0.000110802]; % t2_0 = 0.5 86%

InitialStates = [theta2(1); 0];            % Initial initial states.
Ts            = 0;                                % Time-continuous system.

% 0.074, 0.00012, 0.125, 0.05, 0.04, 0.06, 0.1, 9.81, 4.8, 0.0002

nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
                'Name', 'Rot_Pendulum');



% set(nlgr, 'InputName', 'Voltage', 'InputUnit', 'V',               ...
%           'OutputName', {'Theta 1 Angular position', 'Theta 2 Angular position'}, ...
%           'OutputUnit', {'rad', 'rad'},                         ...
%           'TimeUnit', 's');
% 
% nlgr = setinit(nlgr, 'Name', {'Theta 1 Angular position' 'Theta 1 Angular velocity' ...
%                 'Theta 2 Angular position' 'Theta 2 Angular velocity'});
% nlgr = setinit(nlgr, 'Unit', {'rad' 'rad/s' 'rad' 'rad/s'});
% nlgr = setpar(nlgr, 'Name', {'K_e' 'm1' 'm2' 'l1' 'l2', 'g'});
% nlgr = setpar(nlgr, 'Unit', {'Vs/rad' 'kg' 'kg' 'm' 'm' 'm/s^2'});
nlgr = setpar(nlgr, 'Minimum', {0 -1 0 0.0001});
nlgr = setpar(nlgr, 'Maximum', {0.06 0.1 1 1});

nlgr.SimulationOptions.AbsTol = 1e-6;
nlgr.SimulationOptions.RelTol = 1e-5;
compare(z, nlgr);

nlgr = setinit(nlgr, 'Fixed', {true true}); % Estimate the initial states.
opt = nlgreyestOptions('Display', 'on', 'SearchMethod', 'fmincon');
opt.SearchOptions.MaxIterations = 30;
nlgr = nlgreyest(z, nlgr, opt);

nlgr.Report
fprintf('\n\nThe search termination condition:\n')
nlgr.Report.Termination

% compare(z, nlgr);