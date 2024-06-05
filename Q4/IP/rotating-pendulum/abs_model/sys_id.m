u_struct = load('input_data_chirp_0.25_1_45');
u_data = u_struct.input_data.Data;

theta1_struct = load('theta1_chirp');
theta1_data = theta1_struct.theta1_data.Data;

theta2_struct = load('theta2_chirp');
theta2_data = theta2_struct.theta2_data.Data;

z = iddata([theta1_data, theta2_data], u_data, 0.15, 'Name', 'Rot_Pendulum');
z.InputName = 'Voltage';
z.InputUnit =  'V';
z.OutputName = {'Theta 1 Angular position', 'Theta 2 Angular position'};
z.OutputUnit = {'rad', 'rad'};
z.Tstart = 0;
z.TimeUnit = 's';

% figure('Name', [z.Name ': Voltage input -> Theta 1 Angular position']);
% plot(z(:, 1, 1));   % Plot first input-output pair (Voltage -> Angular position).
% figure('Name', [z.Name ': Voltage input -> Theta 2 Angular position']);
% plot(z(:, 2, 1));   % Plot second input-output pair (Voltage -> Angular velocity).

FileName      = 'rot_pendulum';       % File describing the model structure.
Order         = [2 1 4];           % Model orders [ny nu nx].
%Parameters    = [1; 0.1; 0.04; 0.1; 0.1; 9.81];         % Initial parameters. Np = 6.
Parameters    = [10; 0.1667; 1; 0.981; 41.6667; 0.981];         % Initial parameters. Np = 6.
InitialStates = [0.7090; 0; -0.7085; 0];            % Initial initial states.
Ts            = 0;                 % Time-continuous system.

nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
                'Name', 'Rot_Pendulum');
set(nlgr, 'InputName', 'Voltage', 'InputUnit', 'V',               ...
          'OutputName', {'Theta 1 Angular position', 'Theta 2 Angular position'}, ...
          'OutputUnit', {'rad', 'rad'},                         ...
          'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'Theta 1 Angular position' 'Theta 1 Angular velocity' ...
                'Theta 2 Angular position' 'Theta 2 Angular velocity'});
nlgr = setinit(nlgr, 'Unit', {'rad' 'rad/s' 'rad' 'rad/s'});
nlgr = setpar(nlgr, 'Name', {'K_e' 'm1' 'm2' 'l1' 'l2', 'g'});
nlgr = setpar(nlgr, 'Unit', {'Vs/rad' 'kg' 'kg' 'm' 'm' 'm/s^2'});

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