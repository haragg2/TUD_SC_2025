u_temp = volt.Data(:);
u_data = u_temp(102:2331);

theta2 = detrend(th2.Data(102:2331));
h = volt.Time(2) - volt.Time(1);

z = iddata(theta2, u_data, h, 'Name', 'Rot_Pendulum');
z.InputName = 'Voltage';
z.InputUnit =  'V';
z.OutputName = {'Theta 2 Angular position'};
z.OutputUnit = {'rad'};
z.Tstart = 0;
z.TimeUnit = 's';

figure('Name', [z.Name ': Voltage input -> Theta 2 Angular position']);
plot(z(:));     % Plot input-output pair (Voltage -> Angular velocity).

FileName      = 'rot_pendulum_theta2';       % File describing the model structure.
Order         = [1 1 2];                     % Model orders [ny nu nx].
Parameters    = [0.05914, 0.0563237, 5.65455e-05,  0.000110802];

InitialStates = [theta2(1); 0];            % Initial initial states.
Ts            = 0;                         % Time-continuous system.


nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, 'Name', 'Rot_Pendulum');

set(nlgr, 'InputName', 'Voltage', 'InputUnit', 'V',...
          'OutputName', {'Theta 2 Angular position'},...
          'OutputUnit', {'rad'}, 'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'Theta 2 Angular position' 'Theta 2 Angular velocity'});
nlgr = setinit(nlgr, 'Unit', {'rad' 'rad/s'});
nlgr = setpar(nlgr, 'Name', {'m2' 'c2' 'b2', 'I2'});
nlgr = setpar(nlgr, 'Unit', {'kg' 'm' '' 'kgm^2'});
nlgr = setpar(nlgr, 'Minimum', {0 -0.1 0 0.0001});
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