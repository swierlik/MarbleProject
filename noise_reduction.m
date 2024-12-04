% noise_reduction.m

% Load the data
load('Data/lorenzData.mat') % Contains 'sol', 't', 'dt'
load('Data/systemData.mat') % Contains 'V_x', 'V_y', 'A_x', 'B_x', 'A_y', 'B_y', 'xReg', 'yReg', 'r', 'tspan', 'dt'

load('Data/lorenzDataStochastic.mat') % Contains 'sol2', 't2', 'dt2'
load('Data/systemDataStochastic.mat') % Contains 'V_x2', 'V_y2', 'A_x2', 'B_x2', 'A_y2', 'B_y2', 'xReg2', 'yReg2', 'r2', 'tspan2', 'dt2'

% Extract x, y, z from sol
x_original = sol(:,1);
y_original = sol(:,2);
z_original = sol(:,3);

% Extract x, y, z from sol2
x_noisy = sol2(:,1);
y_noisy = sol2(:,2);
z_noisy = sol2(:,3);

% Define L (exclude initial and final transients)
L = 300:length(xReg)-300;

% ------------------ Noise Reduction on Stochastic Data ------------------

% Universal parameters (from optimization results)
window_size = 11; % Moving Average
order = 2; % Savitzky-Golay Order
framelen = 5; % Savitzky-Golay Frame Length
wd_level = 5; % Wavelet Denoising Level

% Apply noise reduction to x's
x_movmean = movmean(x_noisy, window_size); % Moving Average
x_sg = sgolayfilt(x_noisy, order, framelen); % Savitzky-Golay
x_wd = wdenoise(x_noisy, wd_level); % Wavelet Denoising

% ------------------ Generate HAVOK System Data ------------------
[V_movmean, A_x_movmean, B_x_movmean, xReg_movmean] = getSystem(x_movmean, 100, r, dt, tspan);
[V_sg, A_x_sg, B_x_sg, xReg_sg] = getSystem(x_sg, 100, r, dt, tspan);
[V_wd, A_x_wd, B_x_wd, xReg_wd] = getSystem(x_wd, 100, r, dt, tspan);

% ------------------ Reconstruct and simulate all systems ------------------
L = 1:min(length(tspan), size(xReg, 1));

% Original
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Noisy
sys_x2 = ss(A_x2, B_x2, eye(r-1), 0*B_x2);  % System matrices for x
[y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L, r2), dt2*(L-1), xReg2(1, 1:r2-1));

% Moving Average
sys_x_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);  % System matrices for x
[y_sim_x_movmean, t_sim_x_movmean] = lsim(sys_x_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));

% Savitzky-Golay
sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);  % System matrices for x
[y_sim_x_sg, t_sim_x_sg] = lsim(sys_x_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));

% Wavelet Denoising
sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);  % System matrices for x
[y_sim_x_wd, t_sim_x_wd] = lsim(sys_x_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));

% ------------------ Compute Errors ------------------
error_x = sqrt(sum((V_x(L,1:r-1) - y_sim_x(L,1:r-1)).^2, 2));
error_x_noisy = sqrt(sum((V_x(L,1:r-1) - y_sim_x2(L,1:r-1)).^2, 2));
error_x_movmean = sqrt(sum((V_x(L,1:r-1) - y_sim_x_movmean(L,1:r-1)).^2, 2));
error_x_sg = sqrt(sum((V_x(L,1:r-1) - y_sim_x_sg(L,1:r-1)).^2, 2));
error_x_wd = sqrt(sum((V_x(L,1:r-1) - y_sim_x_wd(L,1:r-1)).^2, 2));

% ------------------ Plot Results ------------------

% Figure 1: Original and noisy forcing vectors
figure;
subplot(2, 2, 1);
plot3(x_original, y_original, z_original, 'b', 'LineWidth', 1);
hold on;
plot3(x_noisy, y_noisy, z_noisy, 'r', 'LineWidth', 1);
title('Original vs Noisy System');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Original', 'Noisy');
grid on;

subplot(2, 2, 2);
plot(x_noisy, 'r', 'LineWidth', 1);
title('Noisy Data');
xlabel('Time');
ylabel('X');
grid on;

subplot(2, 2, 3);
plot(x_movmean, 'g', 'LineWidth', 1);
title('Denoised (Moving Average)');
xlabel('Time');
ylabel('X');
grid on;

subplot(2, 2, 4);
plot(x_sg, 'm', 'LineWidth', 1);
hold on;
plot(x_wd, 'k', 'LineWidth', 1);
title('Denoised (SG & Wavelet)');
xlabel('Time');
ylabel('X');
legend('Savitzky-Golay', 'Wavelet');
grid on;

% Figure 2: Delay Embedded Attractors
figure;
subplot(2, 2, 1);
plot3(V_x(L,1), V_x(L,2), V_x(L,3), 'b', 'LineWidth', 1);
hold on;
plot3(V_x2(L,1), V_x2(L,2), V_x2(L,3), 'r', 'LineWidth', 1);
title('Original vs Noisy Delay Embedded Attractor');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
legend('Original', 'Noisy');
grid on;

subplot(2, 2, 2);
plot3(y_sim_x_movmean(:,1), y_sim_x_movmean(:,2), y_sim_x_movmean(:,3), 'g', 'LineWidth', 1);
title('Reconstructed (Moving Average)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
grid on;

subplot(2, 2, 3);
plot3(y_sim_x_sg(:,1), y_sim_x_sg(:,2), y_sim_x_sg(:,3), 'm', 'LineWidth', 1);
title('Reconstructed (Savitzky-Golay)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
grid on;

subplot(2, 2, 4);
plot3(y_sim_x_wd(:,1), y_sim_x_wd(:,2), y_sim_x_wd(:,3), 'k', 'LineWidth', 1);
title('Reconstructed (Wavelet)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
grid on;

% Figure 3: Error Comparison
figure;
plot(tspan(L), error_x, 'b', 'LineWidth', 1);
hold on;
plot(tspan(L), error_x_noisy, 'r', 'LineWidth', 1);
plot(tspan(L), error_x_movmean, 'g', 'LineWidth', 1);
plot(tspan(L), error_x_sg, 'm', 'LineWidth', 1);
plot(tspan(L), error_x_wd, 'k', 'LineWidth', 1);
title('Error Comparison');
xlabel('Time');
ylabel('Error');
legend('Original','Noisy' , 'Moving Average', 'Savitzky-Golay', 'Wavelet');
grid on;

% ------------------ End of Script ------------------
