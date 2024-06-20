close all;
clear all;

% Set Latex interpreter for plots
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% - Calibration
load("data_graphs\Calibration\Before_Calib.mat");
calib_theta_b = calib_theta(1:900, :);
load("data_graphs\Calibration\After_Calib.mat");
calib_theta_a = calib_theta(1:900, :);
time = 0:0.01:(length(calib_theta_a)-1)*0.01;

figure;
sgtitle({'Calibration of Angular Positions'});
subplot(2, 1, 1);
hold on;
plot(time, calib_theta_b(:, 1), LineWidth=1.2);
plot(time, calib_theta_a(:, 1), LineWidth=1.2);
hold off;
ylabel("$\theta_1$ ($rad$)");
xlabel("Time ($s$)");
grid on;
subplot(2, 1, 2);
hold on;
plot(time, calib_theta_b(:, 2), LineWidth=1.2);
plot(time, calib_theta_a(:, 2), LineWidth=1.2);
hold off;
ylabel("$\theta_2$ ($rad$)");
xlabel("Time ($s$)");
grid on;
legend('Before', 'After', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');

%% Kalman

data = get_data_lqi('data_graphs\lqi\run1_dist.mat', 5000);
figure;
sgtitle("Kalman Filter Performance");
subplot(4, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_offset, 'LineWidth', 1.1);
hplots(2) = plot(data.time, data.th1_hat, 'LineWidth', 1.1);
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(4, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_dot, 'LineWidth', 1.1);
hplots(2) = plot(data.time, data.th1_dot_hat, 'LineWidth', 1.1);
ylabel('$$\dot{\theta_{1}}$$ ($rad/s$)');
hold off;
subplot(4, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_offset, 'LineWidth', 1.1);
hplots(2) = plot(data.time, data.th2_hat, 'LineWidth', 1.1);
ylabel('$\theta_{2}$ ($rad$)');
hold off;
subplot(4, 1, 4);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_dot, 'LineWidth', 1.1);
hplots(2) = plot(data.time, data.th2_dot_hat, 'LineWidth', 1.1);
ylabel('$$\dot{\theta_{2}}$$ ($rad/s$)');
hold off;
xlabel('Time ($s$)');
legend('Actual', 'Estimated', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');

%% - lqr

data = get_data_lqi('data_graphs\lqr\run1_no_dist.mat', 5000);
figure;
sgtitle("Regulation LQR");
subplot(3, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(3, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(3, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
ylabel('$u$ ($V$)');
hold off;
xlabel('Time ($s$)');


%% - lqi

data = get_data_lqi('data_graphs\lqi\run2_no_dist.mat', 5000);
figure;
sgtitle("Regulation LQI");
subplot(4, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(4, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(4, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
ylabel('$u$ ($V$)');
hold off;
K = [0.0187; -0.0187];
int_act = data.u_int1*K(1) + data.u_int2*K(2);
subplot(4, 1, 4);
hold on;
grid on;
hplots(1) = plot(data.time, int_act, 'LineWidth', 1.1);
ylabel('$K_Ix_I$ ($V$)');
hold off;
xlabel('Time ($s$)');

data = get_data_lqi('data_graphs\lqi\run2_dist_clockwise.mat', 5000);
figure;
sgtitle("Regulation LQI with External Disturbance");
subplot(4, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(4, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(4, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
ylabel('$u$ ($V$)');
hold off;
K = [0.0187; -0.0187];
int_act = data.u_int1*K(1) + data.u_int2*K(2);
subplot(4, 1, 4);
hold on;
grid on;
hplots(1) = plot(data.time, int_act, 'LineWidth', 1.1);
ylabel('$K_Ix_I$ ($V$)');
hold off;
xlabel('Time ($s$)');

%% MPC

data = get_data_mpc('data_graphs\Regulation MPC\run1_nointegrator.mat', 5000);
figure;
sgtitle("Regulation MPC without Integrator");
subplot(3, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(3, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(3, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
% yline(mean(data.th1_dot_hat(500:end)), 'r');
ylabel('$u$ ($V$)');
hold off;
xlabel('Time ($s$)');


data = get_data_mpc('data_graphs\Regulation MPC\run1_nodist.mat', 5000);
figure;
sgtitle("Regulation MPC with Integrator");
subplot(3, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(3, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(3, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
% yline(mean(data.th1_dot_hat(500:end)), 'r');
ylabel('$u$ ($V$)');
hold off;
xlabel('Time ($s$)');

data = get_data_mpc('data_graphs\Regulation MPC\run1_dist.mat', 5000);
figure;
sgtitle("Regulation MPC with External Disturbance");
subplot(3, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(3, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(3, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
% yline(mean(data.th1_dot_hat(500:end)), 'r');
ylabel('$u$ ($V$)');
hold off;
xlabel('Time ($s$)');


%% MPC - Reference

data = get_data_mpc('data_graphs\Reference Tracking MPC\Run2_NoDist_pi6.mat', 5000);
figure;
sgtitle("Constant Reference Tracking MPC");
subplot(3, 1, 1);
hold on;
grid on;
hplots(1) = plot(data.time, data.th1_og, 'LineWidth', 1.1);
yline(mean(data.th1_og(500:end)), 'r');
ylabel('$\theta_{1}$ ($rad$)');
hold off;
subplot(3, 1, 2);
hold on;
grid on;
hplots(1) = plot(data.time, data.th2_og, 'LineWidth', 1.1);
yline(mean(data.th2_og(500:end)), 'r');
ylabel('$\theta_{2}$ ($rad$)');
hold off;
legend('Data', 'Mean', 'Location','southoutside', 'Orientation','horizontal', 'Box', 'Off');
subplot(3, 1, 3);
hold on;
grid on;
hplots(1) = plot(data.time, data.u, 'LineWidth', 1.1);
% yline(mean(data.th1_dot_hat(500:end)), 'r');
ylabel('$u$ ($V$)');
hold off;
xlabel('Time ($s$)');


%% Funtions
function data= get_data_lqi(filename, len)
    lqi_data = load(filename);
    data_t = lqi_data.data.extractTimetable;
    data_t = data_t(1:len, :);

    data.time = 0:0.01:(len-1)*0.01;
    data.u = data_t.("Control Input1:1");
    data.th1_og = data_t.th1_orig;
    data.th2_og = data_t.th2_orig;
    data.th1_dot_hat = data_t.th1_dot_hat;
    data.th1_hat = data_t.th1_hat;
    data.th2_dot_hat = data_t.th2_dot_hat;
    data.th2_hat = data_t.th2_hat;
    data.th1_dot = data_t.th1_dot;
    data.th2_dot = data_t.th2_dot;
    data.u_int1 = data_t.("u_int(1)");
    data.u_int2 = data_t.("u_int(2)");
    data.th1_offset = data_t.th1_offset;
    data.th2_offset = data_t.th2_offset;
end

function data= get_data_mpc(filename, len)
    mpc_data = load(filename);
    data_t = mpc_data.data.extractTimetable;
    data_t = data_t(1:len, :);

    data.time = 0:0.01:(len-1)*0.01;
    data.u = data_t.u_sum;
    data.u_mpc = data_t.u_mpc;
    data.ref_err1 = data_t.("ref_error(1)");
    data.ref_err2 = data_t.("ref_error(2)");
    data.th1_og = data_t.th1_orig;
    data.th2_og = data_t.th2_orig;
    data.th1_dot_hat = data_t.th1_dot_hat;
    data.th1_hat = data_t.th1_hat;
    data.th2_dot_hat = data_t.th2_dot_hat;
    data.th2_hat = data_t.th2_hat;
    data.th1_dot = data_t.th1_dot;
    data.th2_dot = data_t.th2_dot;
    data.u_int1 = data_t.("u_int(1)");
    data.u_int2 = data_t.("u_int(2)");
    data.th1_offset = data_t.th1_offset;
    data.th2_offset = data_t.th2_offset;
end