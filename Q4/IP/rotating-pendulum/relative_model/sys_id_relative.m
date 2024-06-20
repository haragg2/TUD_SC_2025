all_data = load('all_0.25_1_02.mat');
u_data = all_data.volt.Data(:);
time = all_data.volt.Time;

theta1 = detrend(all_data.th1.Data(:));
theta2 = detrend(all_data.th2.Data(:));

% u_data = volt.Data(:);
% theta1 = detrend(th1.Data);
% theta2 = detrend(th2.Data);
% h = volt.Time(2) - volt.Time(1);

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

FileName      = 'rot_pendulum_relative';       % File describing the model structure.
Order         = [2 1 4];           % Model orders [ny nu nx].
Parameters    = [0.02, 0.0954944, 0.015, 29.9204, 187.247]; % with limits
InitialStates = [theta1(1); 0; theta2(1); 0];               % Initial initial states.
Ts            = 0;                                          % Time-continuous system.

nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
                'Name', 'Rot_Pendulum');
set(nlgr, 'InputName', 'Voltage', 'InputUnit', 'V',               ...
          'OutputName', {'Theta 1 Angular position', 'Theta 2 Angular position'}, ...
          'OutputUnit', {'rad', 'rad'},                         ...
          'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'Theta 1 Angular position' 'Theta 1 Angular velocity' ...
                'Theta 2 Angular position' 'Theta 2 Angular velocity'});
nlgr = setinit(nlgr, 'Unit', {'rad' 'rad/s' 'rad' 'rad/s'});
nlgr = setpar(nlgr, 'Minimum', {0.02 0.08 1e-3 -50 -100});
nlgr = setpar(nlgr, 'Maximum', {0.3 0.15 0.015 50 200});

nlgr.SimulationOptions.AbsTol = 1e-3;
nlgr.SimulationOptions.RelTol = 1e-3;
compare(z, nlgr);

nlgr = setinit(nlgr, 'Fixed', {false false false false}); % Estimate the initial states.
opt = nlgreyestOptions('Display', 'on');
% opt = nlgreyestOptions('Display', 'on', 'SearchMethod', 'fmincon');
opt.SearchOptions.MaxIterations = 15;
nlgr = nlgreyest(z, nlgr, opt);

nlgr.Report
fprintf('\n\nThe search termination condition:\n')
nlgr.Report.Termination

compare(z, nlgr);