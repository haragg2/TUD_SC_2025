close all;
%clear all;

%load Assignment_Data_SC42145.mat;

%% Part 1.1 to 1.3


% PID controller parameters
Kp = 0.27; 
Td = 1;

% s is the Laplace variable
s = tf('s');

% PID Controller
K_PI_SISO = Kp * (1 + Td*s) / s;

% Assuming FWT is a 3x2 state-space model, extract the appropriate
% transfer function for the beta to omega path. 
G_beta_to_omega = -tf(FWT(1,1));

% Obtain zeros and poles of the transfer function
z_G = zero(G_beta_to_omega);
p_G = pole(G_beta_to_omega);

% Open-loop transfer function
L_beta_to_omega = G_beta_to_omega * K_PI_SISO;

% Closed-loop transfer function
CLTF_beta_to_omega = feedback(L_beta_to_omega, 1);

% Sensitivity function
S_beta_to_omega = 1 - CLTF_beta_to_omega;

% Plotting
figure(1);
bode(G_beta_to_omega, L_beta_to_omega, S_beta_to_omega);
title('Bode Plots of G, L and S');
legend('G_{\beta\rightarrow\omega}', 'L_{\beta\rightarrow\omega}', 'S_{\beta\rightarrow\omega}');
grid on; % Enable grid
print -depsc G_K_S_Bode

figure(2);
rlocus(G_beta_to_omega);
title('Root Locus of G_{\beta\rightarrow\omega}');
print -depsc rlocus

% Step Response with Peak Response and Settling Time
figure(3);
step(CLTF_beta_to_omega);
grid on;
print -depsc KI_Step

% Computing Peak Response and Settling Time
info = stepinfo(CLTF_beta_to_omega);
peakResponse = info.Peak;
settlingTime = info.SettlingTime;

% Displaying Peak Response and Settling Time
disp(['SISO Peak Response: ', num2str(peakResponse)]);
disp(['SISO Settling Time: ', num2str(settlingTime), ' seconds']);

% Check if the closed-loop system is stable
zc = zero(CLTF_beta_to_omega);
pc = pole(CLTF_beta_to_omega);

isStable = all(real(pc) < 0);
if isStable
    disp('The SISO closed-loop system is stable.');
else
    disp('The SISO closed-loop system is unstable.');
end

bw_cl_beta_to_omega = bandwidth(CLTF_beta_to_omega);

%% Part 1.4 disturbance rejection

warning off

% Define the inputs and outputs for FWT and K_PI_SISO
% Assuming FWT has 3 inputs (V, Beta, Tau) and outputs the FWT state
% Assuming K_PI_SISO takes one input and outputs one control action

K_PI_SISO.InputName = 'e';
K_PI_SISO.OutputName = 'u';

S1 = sumblk("e = ref - Omega (rad/s)");
S2 = sumblk("Beta (deg) = -u"); % The plant expects the negative input from the controller
% (because conotroller was designed for the -G plant)


% Connect the systems
inputs = {'ref', 'V (m/s)'};
outputs = {'Omega (rad/s)', 'z (m)', 'Beta (deg)'};
T = connect(FWT, K_PI_SISO, S1, S2, inputs, outputs);

warning on

% Closed-loop system with SISO controller
CL_sisocontroller = minreal(T);

figure();
step(CL_sisocontroller);
grid on;
print -depsc DistRej_Step

%% Part 2

mimo_plant = minreal(balreal(tf(FWT(1:2, 1:2))));

% Compute the RGA
w = 0; %0.3*2*pi;   %bandwidth
s_p = 1j*w;
G_w = evalfr(mimo_plant, s_p);
RGA = G_w.*transpose(pinv(G_w));
RGA_abs = abs(G_w.*pinv(G_w)');

% Compute poles
z_mimo = tzero(mimo_plant);
p_mimo = pole(mimo_plant);
figure();
pzmap(mimo_plant);
title('Pole-zero map of the MIMO plant (G)');
print -depsc MIMO_pole_zero_2.2.eps

% Compute weighting functions
w = 0.3*2*pi;
Ai = 10^-4;
Wp_11 = (s/3 + w)/(s + w*Ai);
Wp = [Wp_11 0; 0 0.2;];
Wu = [0.01 0; 0 (5*10^-3*s^2 + 7*10^-4*s + 5*10^-5)/(s^2 + 14*10^-4*s + 10^-6)];
Wt = 0;

% Generalized Plant
P11 = [zeros(2); Wp];
P12 = [Wu; -Wp * mimo_plant];
P21 = eye(2);
P22 = -mimo_plant;
P = minreal(balreal([P11 P12; P21 P22]));

[K_MIMO,CL,gamma] = hinfsyn(P,2,2);
K_MIMO = minreal(balreal(K_MIMO)); % Why does converting this to tf mess up with the Nyquist plot???

% Gives warning of name conflict
K_MIMO.InputName = {'Omega (rad/s)', 'z (m)'};
K_MIMO.OutputName = {'Beta (deg)', 'tau_e (Nm)'};

S1 = sumblk("e = ref - Omega (rad/s)");
CL_inf_system = stepResponseSimulationMIMO(K_MIMO, FWT, 1000, S1);
print -depsc H_inf_performance.eps

S_infinity_mimo = inv(tf(minreal(balreal(eye(2) + FWT(:,1:2) * K_MIMO))));
T_infinity_mimo = tf(minreal(balreal((FWT(:,1:2) * K_MIMO) * S_infinity_mimo)));

figure();
bodemag(1/Wp_11);
grid on; % Enable grid
title('1/W_{p11} - Performance weighting function on channel 1');
print -depsc wp11_Bode

figure();
Wp.InputName = {'Omega (rad/s)', 'z (m)'};
Wp.OutputName = {'Omega (rad/s)', 'z (m)'};
bodemag(1/Wp);
grid on; % Enable grid
title('1/W_{p} - Performance weighting function');
print -depsc Wp_bode_2.6.eps

figure();
bodemag(1/Wp, S_infinity_mimo);
grid on; % Enable grid
title('1/Wp and Sensitivity of H-infinity');
print -depsc Wp_H_bode.eps
% figure();
% sigma(S_infinity_mimo);

figure();
Wu.InputName = {'Beta (deg)', 'tau_e (Nm)'};
Wu.OutputName = {'Omega (rad/s)', 'z (m)'};
bodemag(1/Wu);
grid on; % Enable grid
title('1/W_{u} - Controller sensitivity weighting function');
print -depsc Wu_bode_2.6.eps

Hinf_controller_sensitivity = K_MIMO * S_infinity_mimo;
figure();
bodemag(1/Wu, Hinf_controller_sensitivity);
grid on; % Enable grid
title('1/Wu and K*S of H-infinity');
print -depsc Wu_KS_infinity_2.7.eps

% Calculate the number of states
K_nstates = size(K_MIMO.A, 1);
G_nstates = size(mimo_plant.A, 1);
Wp_nstates = size(ss(Wp).A, 1);
Wu_nstates = size(ss(Wu).A, 1);
P_nstates = size(P.A, 1);

oneplusL = tf(eye(2) + mimo_plant * K_MIMO);
det_gen_nyq = oneplusL(1,1) * oneplusL(2,2) - oneplusL(1,2) * oneplusL(2,1);
figure();
nyquist(det_gen_nyq);
print -depsc gen_nyquist_1_2.8.eps

%% Part 3.1
Kp = realp('Kp', 1);
Ki = realp('Ki', 1);
Kd = realp('Kd', 1);
Tf = realp('Tf', 1);

Wp_simple = 0.95*(s + 0.04*pi) / (0.016*pi + s);
C_struct = Kp + Ki/s + (Kd*s) / (Tf*s + 1);
% C_struct = Kp + Ki/s;

% SISO hinfstruct
Wp_siso = Wp_simple;
Wp_siso.u = 'y_act';
Wp_siso.y = 'z1';

G_siso = -1 * tf(FWT(1,1));
G_siso.u = 'u';
G_siso.y = 'y_plant';

C_siso = C_struct;
C_siso.u = 'e';
C_siso.y = 'u';

Sum1 = sumblk('e =  y_act');
Sum2 = sumblk('y_act = -y_plant + y_dist');
Siso_Con = connect(G_siso, Wp_siso, C_siso, Sum1, Sum2, {'y_dist'}, {'z1'});

opt = hinfstructOptions('Display', 'final', 'RandomStart', 5);
N_siso = hinfstruct(Siso_Con, opt);

% Extract controller gains:
Kp_opt = N_siso.Blocks.Kp.Value;
Ki_opt = N_siso.Blocks.Ki.Value;
Kd_opt = N_siso.Blocks.Kd.Value ;
Tf_opt = N_siso.Blocks.Tf.Value ;
Kfb_opt = Kp_opt + Ki_opt /s + ( Kd_opt *s)/( Tf_opt *s +1);
%Kfb_opt = Kp_opt + Ki_opt /s;

figure();
bodemag(1/Wp_simple, inv(tf(minreal(balreal(1 + FWT(1,1) * Kfb_opt)))));
grid on;

% Simulate the system
CL_FS_SISO_system = stepResponseSimulationSISO(Kfb_opt, FWT, 300);
print -depsc siso_fwt_pid.eps


%% Part 3.2

Kp_mimo = realp('Kp', [1 0;1 0]);
Ki_mimo = realp('Ki', [1 0;1 0]);
Kd_mimo = realp('Kd', [1 0;1 0]);
Tf_mimo = realp('Tf', 1);

Kp_mimo.Free(1,2) = false;
Kp_mimo.Free(2,2) = false;

Ki_mimo.Free(1,2) = false;
Ki_mimo.Free(2,2) = false;

Kd_mimo.Free(1,2) = false;
Kd_mimo.Free(2,2) = false;

C_struct_mimo = Kp_mimo + Ki_mimo/s + (Kd_mimo*s) / (Tf_mimo*s + 1);

% MIMO hinfstruct
Wp_mimo = Wp;
Wu_mimo = -Wu;

Wp_mimo.u = 'y_act';
Wp_mimo.y = 'z1';

Wu_mimo.u = 'u';
Wu_mimo.y = 'z2';

G_mimo = mimo_plant;
G_mimo.u = 'u';
G_mimo.y = 'y_plant';

C_mimo = C_struct_mimo;
C_mimo.u = 'e';
C_mimo.y = 'u';

Sum1 = sumblk('e = y_act', 2);
Sum2 = sumblk('y_act = -y_plant + y_dist', 2);
Mimo_Con = connect(G_mimo, Wp_mimo, Wu_mimo, C_mimo, Sum1, Sum2, {'y_dist'}, {'z1', 'z2'});

opt = hinfstructOptions('Display', 'final', 'RandomStart', 5);
N_mimo = hinfstruct(Mimo_Con, opt);

% Extract controller gains:
Kp_opt = N_mimo.Blocks.Kp.Value;
Ki_opt = N_mimo.Blocks.Ki.Value ;
Kd_opt = N_mimo.Blocks.Kd.Value ;
Tf_opt = N_mimo.Blocks.Tf.Value ;

% Simulate MIMO using Fixed Structure Controller
Kfb_opt = Kp_opt + Ki_opt /s+( Kd_opt *s)/( Tf_opt *s +1);

S1 = sumblk("e = ref - Omega (rad/s)");
CL_FS_MIMO_system = stepResponseSimulationMIMO(Kfb_opt, FWT, 300, S1);
print -depsc FC_MIMO_performance.eps

S_FS_mimo = inv(tf(minreal(balreal(eye(2) + FWT(:,1:2) * Kfb_opt))));
T_FS_mimo = tf(minreal(balreal((FWT(:,1:2) * Kfb_opt) * S_FS_mimo)));

figure();
bodemag(S_FS_mimo + T_FS_mimo);
grid on;
title('Sensitivity + Complimentary-sensitivity = S + T for Fixed Structure');
print -depsc FC_MIMO_S_T.eps

figure();
bodemag(S_infinity_mimo, S_FS_mimo);
grid on;
title('Sensitivity of H-infinity and FS');
print -depsc FC_H_MIMO_S.eps

figure();
bodemag(T_infinity_mimo, T_FS_mimo);
grid on;
title('Complimentary-sensitivity of H-infinity and FS');
print -depsc FC_H_MIMO_T.eps

%% Function definitions
function CL_system = stepResponseSimulationSISO(controller, plant, stime)
    % Function to simulate the step response of a connected system
    % controller - the controller system
    % stime - step response simulation time
    % plant - the plant system

    controller.InputName = 'e';
    controller.OutputName = 'u';

    % Define the sum blocks
    S1 = sumblk("e = ref - Omega (rad/s)");
    S2 = sumblk("Beta (deg) = -u");

    % Connect the systems
    inputs = {'ref', 'V (m/s)'};
    outputs = {'Omega (rad/s)', 'z (m)', 'Beta (deg)'};
    T = connect(plant, controller, S1, S2, inputs, outputs);

    % Simplify the closed-loop system
    CL_system = minreal(T);

    % Perform the step response simulation
    figure();
    step(CL_system, stime);
    title('Step Response of the SISO-Loop System');
    grid on;
    xlabel('Time');
    ylabel('Response');

end

function CL_system = stepResponseSimulationMIMO(controller, plant, stime, varargin)
    % Function to simulate the step response with variable number of sum blocks
    % controller - the controller system
    % plant - the plant system
    % stime - step response simulation time
    % varargin - variable number of sum blocks

    controller.Inputname = {'e', 'z (m)'};
    controller.Outputname = {'Beta (deg)', 'tau_e (Nm)'};

    % Define inputs and outputs for the system
    inputs = {'ref', 'V (m/s)'};
    outputs = {'Omega (rad/s)', 'z (m)', 'Beta (deg)', 'tau_e (Nm)'};

    % Define the system connections
    sumBlocks = varargin;  % Extract sum blocks from varargin

    % Connect the systems with the sum blocks
    connectedBlocks = [{plant, controller}, sumBlocks];
    T = connect(connectedBlocks{:}, inputs, outputs);

    % Simplify the closed-loop system
    CL_system = minreal(T);

    % Perform the step response simulation
    figure();
    step(CL_system, stime);
    title('Step Response of the MIMO-Loop System');
    xlabel('Time');
    ylabel('Response');
    grid on;
    
end
