% noiseAccuracyGraph.m
% This program analyzes how different noise levels affect the reconstruction accuracy
% of Lorenz system data using different noise reduction methods

% Clear workspace and close figures
clear all;
close all;

% Load the original deterministic Lorenz data
load('Data/lorenzData.mat');

% Define noise levels to test (logarithmic scale)
noise_levels = logspace(-10, 0, 20); % From 10^-10 to 10^0 (1)

% Initialize arrays to store Mean Square Error for each method
mse_noisy = zeros(size(noise_levels));
mse_movmean = zeros(size(noise_levels));
mse_sg = zeros(size(noise_levels));
mse_wd = zeros(size(noise_levels));


% Define parameters for noise reduction methods
window_size = 7;    % Moving Average
sg_order = 2;       % Savitzky-Golay Order
sg_framelen = 23;   % Savitzky-Golay Frame Length
wd_level = 5;       % Wavelet Denoising Level
wd_wavelet = 'sym4'; % Wavelet Denoising Wavelet

% Extract original x data
x_original = sol(:,1);

% Loop through each noise level
for i = 1:length(noise_levels)
    eta = noise_levels(i);
    fprintf('Processing noise level: %.10g\n', eta);
    
    % Generate noisy data
    sol_noisy = sol;
    sol_noisy(:,1) = sol(:,1) + eta * randn(size(sol,1), 1);
    
    % Extract noisy x data
    x_noisy = sol_noisy(:,1);
    
    % Apply noise reduction methods
    x_movmean = movmean(x_noisy, window_size);
    x_sg = sgolayfilt(x_noisy, sg_order, sg_framelen);
    x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet);
    
    % Generate HAVOK systems for each method
    r = 15; % Rank for HAVOK model
    
    % Process each dataset
    [V_orig, A_x_orig, B_x_orig, xReg_orig, U_orig, E_orig] = getSystem(x_original, 100, r, dt, t);
    [V_noisy, A_x_noisy, B_x_noisy, xReg_noisy, U_noisy, E_noisy] = getSystem(x_noisy, 100, r, dt, t);
    [V_movmean, A_x_movmean, B_x_movmean, xReg_movmean, U_movmean, E_movmean] = getSystem(x_movmean, 100, r, dt, t);
    [V_sg, A_x_sg, B_x_sg, xReg_sg, U_sg, E_sg] = getSystem(x_sg, 100, r, dt, t);
    [V_wd, A_x_wd, B_x_wd, xReg_wd, U_wd, E_wd] = getSystem(x_wd, 100, r, dt, t);
    
    % Define common L range
    L = 300:min([size(xReg_orig,1), size(xReg_noisy,1), size(xReg_movmean,1), size(xReg_sg,1), size(xReg_wd,1)]) - 300;
    
    % Create state space systems
    sys_orig = ss(A_x_orig, B_x_orig, eye(r-1), 0*B_x_orig);
    sys_noisy = ss(A_x_noisy, B_x_noisy, eye(r-1), 0*B_x_noisy);
    sys_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);
    sys_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);
    sys_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);
    
    % Simulate systems
    [y_sim_orig, t_sim_orig] = lsim(sys_orig, xReg_orig(L, r), dt*(L-1), xReg_orig(1, 1:r-1));
    [y_sim_noisy, t_sim_noisy] = lsim(sys_noisy, xReg_noisy(L, r), dt*(L-1), xReg_noisy(1, 1:r-1));
    [y_sim_movmean, t_sim_movmean] = lsim(sys_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));
    [y_sim_sg, t_sim_sg] = lsim(sys_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));
    [y_sim_wd, t_sim_wd] = lsim(sys_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));
    
    % Ensure correct dimensions for reconstruction
    U_orig = U_orig(:, 1:r-1);
    U_noisy = U_noisy(:, 1:r-1);
    U_movmean = U_movmean(:, 1:r-1);
    U_sg = U_sg(:, 1:r-1);
    U_wd = U_wd(:, 1:r-1);
    
    E_orig = E_orig(1:r-1, 1:r-1);
    E_noisy = E_noisy(1:r-1, 1:r-1);
    E_movmean = E_movmean(1:r-1, 1:r-1);
    E_sg = E_sg(1:r-1, 1:r-1);
    E_wd = E_wd(1:r-1, 1:r-1);
    
    % Reconstruct signals
    x_reconstructed_orig = U_orig * E_orig * y_sim_orig.';
    x_reconstructed_noisy = U_noisy * E_noisy * y_sim_noisy.';
    x_reconstructed_movmean = U_movmean * E_movmean * y_sim_movmean.';
    x_reconstructed_sg = U_sg * E_sg * y_sim_sg.';
    x_reconstructed_wd = U_wd * E_wd * y_sim_wd.';
    
    % Calculate MSE - Use original x as ground truth
    ref_x = x_original(L);
    
    % Cut to ensure matching lengths (if necessary)
    min_length = min([length(ref_x), size(x_reconstructed_orig, 2), size(x_reconstructed_noisy, 2), ...
                     size(x_reconstructed_movmean, 2), size(x_reconstructed_sg, 2), ...
                     size(x_reconstructed_wd, 2)]);
    
    ref_x = ref_x(1:min_length)';
    x_reconstructed_orig = x_reconstructed_orig(1, 1:min_length);
    x_reconstructed_noisy = x_reconstructed_noisy(1, 1:min_length);
    x_reconstructed_movmean = x_reconstructed_movmean(1, 1:min_length);
    x_reconstructed_sg = x_reconstructed_sg(1, 1:min_length);
    x_reconstructed_wd = x_reconstructed_wd(1, 1:min_length);
    
    % Calculate MSE
    mse_noisy(i) = mean((ref_x - x_reconstructed_noisy).^2);
    mse_movmean(i) = mean((ref_x - x_reconstructed_movmean).^2);
    mse_sg(i) = mean((ref_x - x_reconstructed_sg).^2);
    mse_wd(i) = mean((ref_x - x_reconstructed_wd).^2);
end

% Create log-log plot of noise level vs. MSE
figure;
loglog(noise_levels, mse_noisy, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(noise_levels, mse_movmean, 'gs-', 'LineWidth', 2, 'MarkerSize', 8);
loglog(noise_levels, mse_sg, 'bd-', 'LineWidth', 2, 'MarkerSize', 8);
loglog(noise_levels, mse_wd, 'kx-', 'LineWidth', 2, 'MarkerSize', 8);

% Add reference lines for noise^1 and noise^2 scaling
ref_x = logspace(-10, 0, 100);
ref_y1 = 10 * ref_x;  % Linear scaling (noise^1)
ref_y2 = 10 * ref_x.^2;  % Quadratic scaling (noise^2)
loglog(ref_x, ref_y1, 'k--', 'LineWidth', 1);
loglog(ref_x, ref_y2, 'k:', 'LineWidth', 1);

% Add labels
xlabel('Noise Level (η)', 'FontSize', 14);
ylabel('Mean Square Error', 'FontSize', 14);
title('Reconstruction Error vs. Noise Level', 'FontSize', 16);
legend('No Filtering', 'Moving Average', 'Savitzky-Golay', 'Wavelet', ...
       'Linear Scaling (∝ η)', 'Quadratic Scaling (∝ η²)', ...
       'Location', 'NorthWest');
grid on;
ylim([1e-10 1e1]); % Set y-axis limits from 10^-10 to 10^1

% Optional: Create a second plot showing the relative improvement
figure;
semilogx(noise_levels, mse_noisy ./ mse_noisy, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogx(noise_levels, mse_movmean ./ mse_noisy, 'gs-', 'LineWidth', 2, 'MarkerSize', 8);
semilogx(noise_levels, mse_sg ./ mse_noisy, 'bd-', 'LineWidth', 2, 'MarkerSize', 8);
semilogx(noise_levels, mse_wd ./ mse_noisy, 'kx-', 'LineWidth', 2, 'MarkerSize', 8);

% Add reference line for no improvement
ref_line = ones(size(noise_levels));
semilogx(noise_levels, ref_line, 'k--', 'LineWidth', 1);

% Add labels
xlabel('Noise Level (η)', 'FontSize', 14);
ylabel('Relative Error (compared to no filtering)', 'FontSize', 14);
title('Relative Improvement of Noise Reduction Methods', 'FontSize', 16);
legend('No Filtering (baseline)', 'Moving Average', 'Savitzky-Golay', 'Wavelet', ...
       'No Improvement Line', 'Location', 'Best');
grid on;
ylim([1e-10 1e1]); % Set y-axis limits from -10 to 10

% % Save the figures
% saveas(figure(1), 'noise_accuracy_mse.png');
% saveas(figure(2), 'noise_accuracy_relative.png');

% % Save the results to a MAT file for later analysis
% save('noise_accuracy_results.mat', 'noise_levels', 'mse_noisy', 'mse_movmean', 'mse_sg', 'mse_wd');

fprintf('Analysis complete. Results saved to noise_accuracy_results.mat\n');