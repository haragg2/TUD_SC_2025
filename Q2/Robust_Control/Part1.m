close all;
%clear all;

%load Assignment_Data_SC42145.mat;
s = tf('s');
Kp = 0.2552;
%Kp = 1.5;
Td = 4.1;

K_PI = Kp*(1+Td*s)/s;

TFs=tf(FWT);
[num_TF11,den_TF11]=tfdata(-TFs(1,1),'v');
TF11 = -TFs(1,1);
z = zero(TF11);
p = pole(TF11);

OLTF = K_PI*TF11;
CLTF = feedback(OLTF,1);
S = 1 - CLTF;
figure(1);
%bode(TF11,K*TF11,CLTF,S);
bode(TF11,OLTF,S);
figure(2);
rlocus(TF11);

%sisotool(TF11);

figure(3);
step(CLTF);
bw = bandwidth(CLTF);
zc = zero(CLTF);
pc = pole(CLTF);
% bw_ol = bandwidth(TF11)
% 
% sys_dist = TFs(1,3);
% CLTF_dist = feedback(K*sys_dist,1);
% S_dist = 1- CLTF_dist;
% figure(4);
% bode(CLTF_dist,S_dist);
% figure(5);
% step(CLTF_dist,1000);

mimo_plant = -minreal(balreal(TFs(1:2, 1:2)));
w = 0.3*2*pi;
Ai = 10^-4;
s_p = j*w;
G_w = evalfr(mimo_plant,s_p);
RGA = G_w.*transpose(pinv(G_w));
RGA_abs = abs(G_w.*pinv(G_w)');

z_mimo = tzero(mimo_plant);
p_mimo = pole(mimo_plant);
figure(6);
pzmap(mimo_plant);

Wp_11 = (s/3 + w)/(s + w*Ai);
Wp = [Wp_11 0; 0 0.2;];
Wu = [0.01 0; 0 (5*10^-3*s^2 + 7*10^-4*s + 5*10^-5)/(s^2 + 14*10^-4*s + 10^-6)];
Wt=0;
figure();
bode(Wu);

[K3,CL3,GAM3,INFO3]=mixsyn(mimo_plant,Wp,Wu,Wt);

P11 = [Wp; zeros(2);];
P12 = [Wp*mimo_plant; Wu;];
P21 = -eye(2);
P22 = -mimo_plant;
P = [P11 P12; P21 P22];
[K,CL,gamma] = hinfsyn(P,2,2) 

sisocontroller=Kp*tf([Td 1], [1 0]); % simple controller with positive feedback

TF_V_w = TFs(1,3);
TF_V_z = TFs(2,3);

systemnames ='mimo_plant Wp Wu';% TF_V_w TF_V_z'; % Define systems
inputvar ='[w(2); u(2)]'; % Input generalized plant
%input_to_TF_V_w= '[u]';
%input_to_TF_V_z= '[u]';
input_to_mimo_plant= '[u]';
input_to_Wu= '[u]';
input_to_Wp= '[w+mimo_plant]';
outputvar= '[Wp; Wu; -mimo_plant-w]'; % Output generalized plant
sysoutname='P1';
sysic;
[K2,CL2,GAM2,INFO2] = hinfsyn(P1,2,2); % Hinf design

oneplusL = tf(eye(2) + K2*mimo_plant);
det_gen_nyq = oneplusL(1,1) * oneplusL(1,2) - oneplusL(2,1) * oneplusL(2,2);
figure();
nyquist(det_gen_nyq);

mimo_CLTF = feedback(mimo_plant*K3,eye(2));
figure();
step(mimo_CLTF(2,2));
figure();
step(mimo_CLTF(2,1));
figure();
step(mimo_CLTF);

 %% disturbance rejection

warning off

systemnames ='FWT';     % The full wind turbine model

inputvar ='[V; Beta; Tau]';    % Input generalized plant

input_to_FWT= '[Beta; Tau; V]';

outputvar= '[FWT; Beta; Tau; FWT]';    % Output generalized plant also includes the inputs

sysoutname='Gsim';
cleansysic = 'yes';

sysic;

warning on

% CL_sisocontroller=minreal(lft(Gsim(1:end-1,1:end-1),sisocontroller)); % SISO controller
CL_mimocontroller=minreal(lft(Gsim, K)) % MIMO controller

figure();

%step(CL_sisocontroller); % simple code for a step on the wind
step(CL_mimocontroller);