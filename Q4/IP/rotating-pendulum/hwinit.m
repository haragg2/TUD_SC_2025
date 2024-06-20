%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSCS FPGA interface board: init and I/O conversions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gains and offsets
daoutoffs = [0.00];                   % output offset
daoutgain = 1*[-6];                   % output gain

% Sensor calibration:
adinoffs = -[3.7653 1.213];
adingain = -[1.207 1.207];

adinoffs = [adinoffs 0 0 0 0 0];    % input offset
adingain = [adingain 1 1 1 1 1];     % input gain (to radians)

h = 0.01;

%% generate signals
t = 0:h:60;
chirp_out = chirp(t,.25,t(end),1);
chirp_out = 0.2*timeseries(chirp_out,t);
