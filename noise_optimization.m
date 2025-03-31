% denoise_parameter_optimization.m

% Load the data
load('Data/lorenzData.mat') % Contains 'sol', 't', 'dt'

% Extract x from sol
x_original = sol(:,1);

% Parameter ranges to test
variance = 1; % Fixed noise level for testing
tspan = dt:dt:50;
r = 5; % HAVOK rank parameter

% Parameter ranges
window_sizes = [3, 5, 7, 9, 11, 13, 15];
sg_orders = [2,3,4];
sg_framelens = [11:4:111]; % Must be odd and > order
wd_levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
wd_wavelets = {'haar', 'db4', 'sym4', 'coif3', 'bior3.5', 'dmey'};

% Add Gaussian Noise to x_original
x_noisy = x_original + variance*randn(size(x_original));

% Initialize results containers
results = struct();
results.moving_avg = struct('window_size', [], 'error', Inf);
results.sg = struct('order', [], 'framelen', [], 'error', Inf);
results.wavelet = struct('level', [], 'wavelet', '', 'error', Inf);

% Function to compute reconstruction error
function error_val = compute_reconstruction_error(x_denoised, x_original, r, dt, tspan)
    % Get HAVOK system for the denoised data
    [V, A, B, xReg, U, E] = getSystem(x_denoised, 100, r, dt, tspan);
    
    % Simulate the system
    L = 1:min(length(tspan), size(xReg, 1));
    sys = ss(A, B, eye(r-1), 0*B);
    [y_sim, t_sim] = lsim(sys, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));
    
    % Ensure U and E match Koopman mode rank
    U = U(:, 1:r-1);
    E = E(1:r-1, 1:r-1);
    
    % Reconstruct
    x_reconstructed = U * E * y_sim.';
    x_reconstructed = dehankelize(x_reconstructed);
    
    % Compute error relative to original data
    error_val = mean((x_original(L) - x_reconstructed(L)).^2);
end

% -------- Optimize Moving Average Parameters --------
fprintf('Optimizing Moving Average parameters...\n');
best_ma_error = Inf;
best_window = 0;

for window_size = window_sizes
    fprintf('Testing window size: %d\n', window_size);
    
    % Apply Moving Average
    x_movmean = movmean(x_noisy, window_size);
    
    % Compute error
    error_val = compute_reconstruction_error(x_movmean, x_original, r, dt, tspan);
    
    if error_val < best_ma_error
        best_ma_error = error_val;
        best_window = window_size;
        fprintf('New best window size: %d (Error: %.8f)\n', window_size, error_val);
    end
end

results.moving_avg.window_size = best_window;
results.moving_avg.error = best_ma_error;

% -------- Optimize Savitzky-Golay Parameters --------
fprintf('\nOptimizing Savitzky-Golay parameters...\n');
best_sg_error = Inf;
best_order = 0;
best_framelen = 0;

for order = sg_orders
    for framelen = sg_framelens
        % Skip invalid combinations (framelen must be > order and odd)
        if framelen <= order || mod(framelen, 2) == 0
            continue;
        end
        
        fprintf('Testing order: %d, framelen: %d\n', order, framelen);
        
        % Apply Savitzky-Golay
        x_sg = sgolayfilt(x_noisy, order, framelen);
        
        % Compute error
        error_val = compute_reconstruction_error(x_sg, x_original, r, dt, tspan);
        
        if error_val < best_sg_error
            best_sg_error = error_val;
            best_order = order;
            best_framelen = framelen;
            fprintf('New best SG parameters: order=%d, framelen=%d (Error: %.8f)\n', order, framelen, error_val);
        end
    end
end

results.sg.order = best_order;
results.sg.framelen = best_framelen;
results.sg.error = best_sg_error;

% -------- Optimize Wavelet Denoising Parameters --------
fprintf('\nOptimizing Wavelet Denoising parameters...\n');
best_wd_error = Inf;
best_level = 0;
best_wavelet = '';

for level = wd_levels
    for w = 1:length(wd_wavelets)
        wavelet = wd_wavelets{w};
        fprintf('Testing level: %d, wavelet: %s\n', level, wavelet);
        
        % Apply Wavelet Denoising
        try
            x_wd = wdenoise(x_noisy, level, 'Wavelet', wavelet);
            
            % Compute error
            error_val = compute_reconstruction_error(x_wd, x_original, r, dt, tspan);
            
            if error_val < best_wd_error
                best_wd_error = error_val;
                best_level = level;
                best_wavelet = wavelet;
                fprintf('New best WD parameters: level=%d, wavelet=%s (Error: %.8f)\n', level, wavelet, error_val);
            end
        catch e
            fprintf('Error with wavelet %s at level %d: %s\n', wavelet, level, e.message);
        end
    end
end

results.wavelet.level = best_level;
results.wavelet.wavelet = best_wavelet;
results.wavelet.error = best_wd_error;

% -------- Display Results --------
fprintf('\n\n===== OPTIMIZATION RESULTS =====\n');
fprintf('For noise level: %.4f\n\n', variance);

fprintf('Best Moving Average Parameters:\n');
fprintf('  Window Size: %d\n', results.moving_avg.window_size);
fprintf('  Reconstruction Error: %.8f\n\n', results.moving_avg.error);

fprintf('Best Savitzky-Golay Parameters:\n');
fprintf('  Order: %d\n', results.sg.order);
fprintf('  Frame Length: %d\n', results.sg.framelen);
fprintf('  Reconstruction Error: %.8f\n\n', results.sg.error);

fprintf('Best Wavelet Denoising Parameters:\n');
fprintf('  Decomposition Level: %d\n', results.wavelet.level);
fprintf('  Wavelet: %s\n', results.wavelet.wavelet);
fprintf('  Reconstruction Error: %.8f\n', results.wavelet.error);

% -------- Create Final Reconstructions with Best Parameters --------
% Apply best denoising methods
x_best_movmean = movmean(x_noisy, results.moving_avg.window_size);
x_best_sg = sgolayfilt(x_noisy, results.sg.order, results.sg.framelen);
x_best_wd = wdenoise(x_noisy, results.wavelet.level, 'Wavelet', results.wavelet.wavelet);

% Generate HAVOK systems
[V_x, A_x, B_x, xReg, U_x, E_x] = getSystem(x_original, 100, r, dt, tspan);
[V_x2, A_x2, B_x2, xReg2, U_x2, E_x2] = getSystem(x_noisy, 100, r, dt, tspan);
[V_movmean, A_x_movmean, B_x_movmean, xReg_movmean, U_movmean, E_movmean] = getSystem(x_best_movmean, 100, r, dt, tspan);
[V_sg, A_x_sg, B_x_sg, xReg_sg, U_sg, E_sg] = getSystem(x_best_sg, 100, r, dt, tspan);
[V_wd, A_x_wd, B_x_wd, xReg_wd, U_wd, E_wd] = getSystem(x_best_wd, 100, r, dt, tspan);

% Reconstruct and simulate all systems
L = 1:min(length(tspan), size(xReg, 1));

% Original
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Noisy
sys_x2 = ss(A_x2, B_x2, eye(r-1), 0*B_x2);
[y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L, r), dt*(L-1), xReg2(1, 1:r-1));

% Moving Average
sys_x_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);
[y_sim_x_movmean, t_sim_x_movmean] = lsim(sys_x_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));

% Savitzky-Golay
sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);
[y_sim_x_sg, t_sim_x_sg] = lsim(sys_x_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));

% Wavelet Denoising
sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);
[y_sim_x_wd, t_sim_x_wd] = lsim(sys_x_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));

% Ensure Koopman modes match
U_x = U_x(:, 1:r-1);
U_x2 = U_x2(:, 1:r-1);
U_movmean = U_movmean(:, 1:r-1);
U_sg = U_sg(:, 1:r-1);
U_wd = U_wd(:, 1:r-1);

% Ensure E matches Koopman mode rank
E_x = E_x(1:r-1, 1:r-1);
E_x2 = E_x2(1:r-1, 1:r-1);
E_movmean = E_movmean(1:r-1, 1:r-1);
E_sg = E_sg(1:r-1, 1:r-1);
E_wd = E_wd(1:r-1, 1:r-1);

% Project simulated results back
x_reconstructed = U_x * E_x * y_sim_x.';
x_reconstructed_noisy = U_x2 * E_x2 * y_sim_x2.';
x_reconstructed_movmean = U_movmean * E_movmean * y_sim_x_movmean.';
x_reconstructed_sg = U_sg * E_sg * y_sim_x_sg.';
x_reconstructed_wd = U_wd * E_wd * y_sim_x_wd.';

x_reconstructed = dehankelize(x_reconstructed);
x_reconstructed_noisy = dehankelize(x_reconstructed_noisy);
x_reconstructed_movmean = dehankelize(x_reconstructed_movmean);
x_reconstructed_sg = dehankelize(x_reconstructed_sg);
x_reconstructed_wd = dehankelize(x_reconstructed_wd);

% Error calculation
error_x = (x_original(L) - x_reconstructed(L)).^2;
error_x_noisy = (x_original(L) - x_reconstructed_noisy(L)).^2;
error_x_movmean = (x_original(L) - x_reconstructed_movmean(L)).^2;
error_x_sg = (x_original(L) - x_reconstructed_sg(L)).^2;
error_x_wd = (x_original(L) - x_reconstructed_wd(L)).^2;

% -------- Plot Results --------
% Figure 1: Reconstruction comparison
figure;
hold on;
plot(t_sim_x, x_reconstructed(L), 'b', 'LineWidth', 1.2, 'DisplayName', 'HAVOK Original');
plot(t_sim_x2, x_reconstructed_noisy(L), 'r', 'LineWidth', 1.2, 'DisplayName', 'Noisy HAVOK');
plot(t_sim_x_movmean, x_reconstructed_movmean(L), 'g', 'LineWidth', 1.2, 'DisplayName', 'Moving Avg');
plot(t_sim_x_sg, x_reconstructed_sg(L), 'm', 'LineWidth', 1.2, 'DisplayName', 'Savitzky-Golay');
plot(t_sim_x_wd, x_reconstructed_wd(L), 'c', 'LineWidth', 1.2, 'DisplayName', 'Wavelet');
xlabel('Time');
ylabel('x(t)');
title('Reconstruction of x(t) using HAVOK with Optimized Parameters');
legend('Location', 'best');
grid on;

% Figure 2: Error comparison
figure;
% Calculate average error in sliding windows to smooth the plot
window = 100;
error_x_smooth = movmean(error_x, window);
error_x_noisy_smooth = movmean(error_x_noisy, window);
error_x_movmean_smooth = movmean(error_x_movmean, window);
error_x_sg_smooth = movmean(error_x_sg, window);
error_x_wd_smooth = movmean(error_x_wd, window);

hold on;
plot(tspan(L), error_x_smooth, 'b', 'LineWidth', 1.2, 'DisplayName', 'Original');
plot(tspan(L), error_x_noisy_smooth, 'r', 'LineWidth', 1.2, 'DisplayName', 'Noisy');
plot(tspan(L), error_x_movmean_smooth, 'g', 'LineWidth', 1.2, 'DisplayName', 'Moving Avg');
plot(tspan(L), error_x_sg_smooth, 'm', 'LineWidth', 1.2, 'DisplayName', 'Savitzky-Golay');
plot(tspan(L), error_x_wd_smooth, 'c', 'LineWidth', 1.2, 'DisplayName', 'Wavelet');
xlabel('Time');
ylabel('Squared Error (Smoothed)');
title('Error Comparison with Optimized Parameters');
legend('Location', 'best');
grid on;

% Figure 3: Bar chart of mean errors
figure;
method_names = {'Original', 'Noisy', 'Moving Avg', 'S-G', 'Wavelet'};
mean_errors = [mean(error_x), mean(error_x_noisy), mean(error_x_movmean), mean(error_x_sg), mean(error_x_wd)];
bar(mean_errors);
set(gca, 'XTickLabel', method_names);
title('Mean Squared Error by Method');
ylabel('Mean Squared Error');
grid on;

% Figure 4: Parameter sensitivity analysis
% This would require running the optimization again with different configurations
% For simplicity, we'll just plot the error values from our optimization

% Create a table with all results
fprintf('\nSaving results to optimization_results.mat\n');
save('optimization_results.mat', 'results', 'x_best_movmean', 'x_best_sg', 'x_best_wd', ...
    'x_reconstructed', 'x_reconstructed_noisy', 'x_reconstructed_movmean', 'x_reconstructed_sg', 'x_reconstructed_wd', ...
    'error_x', 'error_x_noisy', 'error_x_movmean', 'error_x_sg', 'error_x_wd');

fprintf('\nDone! Optimized parameters have been found and results saved.\n');