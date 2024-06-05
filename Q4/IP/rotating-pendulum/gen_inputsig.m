
t = 0:h:60;
chirp_out = chirp(t,.5,t(end),2);
chirp_out = 0.3*timeseries(chirp_out,t);
figure
plot(chirp_out)

% sig = zeros(1, size(t,2));
% sig(500:499+0.5/h) = 0.5*ones(1,0.5/h);
% sig = 0*timeseries(sig,t);
% plot(sig)