%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSCS FPGA interface board: init and I/O conversions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gains and offsets
daoutoffs = [0.00];                   % output offset
daoutgain = 1*[-6];                   % output gain

% Sensor calibration:
adinoffs = -[0 0];
adingain = [1 1];

% old calib
% adinoffs = -[3.786 1.243];
% adingain = -[1.207 1.207];
% % % % adinoffs = -[3.7547 1.199];
%latest calib
adinoffs = -[3.7653 1.213];
adingain = -[1.207 1.207];
% daoutoffs = [-0.00002]; 

adinoffs = [adinoffs 0 0 0 0 0];    % input offset
adingain = [adingain 1 1 1 1 1];     % input gain (to radians)

h = 0.01;

%% generate signals
t = 0:h:60;
chirp_out = chirp(t,.25,t(end),1);
chirp_out = 0.2*timeseries(chirp_out,t);
%figure
% plot(chirp_out)

sig = zeros(1, size(t,2));
sig(500:499+0.5/h) = 0.5*ones(1,0.5/h);
sig = 0*timeseries(sig,t);
% plot(sig)