% noise_reduction_modified.m
clear; clc; close all;

% -------- USER DEFINED PARAMETERS --------
TARGET_R = 7;          % Set the desired HAVOK rank for this run
TARGET_VARIANCE = 3; % Set the desired noise level variance for this run
results_filename = 'Data/optimal_denoise_params.mat'; % File with precomputed params
run_optimization_if_missing = true; % Set to true to run optimizer if params not found
% -----------------------------------------

fprintf('Running analysis for r = %d, variance = %g\n', TARGET_R, TARGET_VARIANCE);

% Load the original data
try
    load('Data/lorenzData.mat', 'sol', 't', 'dt'); % Contains 'sol', 't', 'dt'
    x_original = sol(:,1);
    tspan = dt:dt:t(end);
    fprintf('Loaded Lorenz data.\n');
catch ME
    error('Failed to load Data/lorenzData.mat: %s', ME.message);
end

% Add Gaussian Noise to x_original
rng('default'); % Reset RNG for consistent noise generation if script is rerun
rng(1); % Or use a specific seed if needed across runs
x_noisy = x_original + sqrt(TARGET_VARIANCE) * randn(size(x_original));
fprintf('Added Gaussian noise with variance %g.\n', TARGET_VARIANCE);

% --- Parameter Loading ---
params_found = false;
bestParams = struct(); % Initialize empty struct

% Try loading precomputed parameters
if exist(results_filename, 'file')
    fprintf('Loading precomputed parameters from %s...\n', results_filename);
    load(results_filename, 'optimization_results');

    % Construct field names to look for
    r_field = sprintf('r%d', TARGET_R);
    var_field = sprintf('var%g', TARGET_VARIANCE);
    var_field = strrep(var_field, '.', 'p');
    var_field = strrep(var_field, '-', 'neg');

    % Check if the specific combination exists
    if isfield(optimization_results, r_field) && isfield(optimization_results.(r_field), var_field)
        bestParams = optimization_results.(r_field).(var_field);
        params_found = true;
        fprintf('Found precomputed parameters for r=%d, var=%g.\n', TARGET_R, TARGET_VARIANCE);
    else
        fprintf('Precomputed parameters NOT found for r=%d, var=%g.\n', TARGET_R, TARGET_VARIANCE);
    end
else
    fprintf('Precomputed parameter file %s not found.\n', results_filename);
end

% --- Run Optimization if Parameters Not Found ---
if ~params_found
    if run_optimization_if_missing
        fprintf('Running optimization on-the-fly...\n');

        % Define parameter ranges FOR ON-THE-FLY optimization (can be same or different)
        paramRanges.window_sizes = [3, 5, 7, 9, 11, 13, 15];
        paramRanges.sg_orders    = [2, 3, 4];
        paramRanges.sg_framelens = [11:4:111];
        paramRanges.wd_levels    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        paramRanges.wd_wavelets  = {'haar', 'db4', 'sym4', 'coif3', 'bior3.5', 'dmey'};

        % Ensure findOptimalDenoiseParams is available
        if exist('findOptimalDenoiseParams', 'file') ~= 2
             error('findOptimalDenoiseParams.m function not found in the MATLAB path.');
        end

         try
            bestParams = findOptimalDenoiseParams(x_original, x_noisy, TARGET_R, dt, tspan, paramRanges);
            params_found = true; % Parameters are now available
            fprintf('On-the-fly optimization complete.\n');

            % --- Optional: Save newly computed parameters ---
            save_new = input('Save these newly computed parameters? (y/n): ', 's');
            if lower(save_new) == 'y'
                if ~exist('optimization_results','var') % If file didn't exist before
                    optimization_results = struct();
                end
                r_field = sprintf('r%d', TARGET_R);
                var_field = sprintf('var%g', TARGET_VARIANCE);
                var_field = strrep(var_field, '.', 'p');
                var_field = strrep(var_field, '-', 'neg');
                optimization_results.(r_field).(var_field) = bestParams;
                save(results_filename, 'optimization_results');
                fprintf('Saved new parameters to %s\n', results_filename);
            end
            % --- End Optional Save ---

         catch ME
              error('Error during on-the-fly optimization: %s\nCheck findOptimalDenoiseParams.m and its dependencies.', ME.message);
         end

    else
        error('Required parameters not found in %s and on-the-fly optimization is disabled.', results_filename);
    end
end

% --- Assign Parameters ---
if ~params_found || isempty(fieldnames(bestParams))
   error('Failed to obtain denoising parameters.');
end

window_size = bestParams.ma_window_size;
order       = bestParams.sg_order;
framelen    = bestParams.sg_framelen;
wd_level    = bestParams.wd_level;
wd_wavelet  = bestParams.wd_wavelet;

fprintf('\nUsing Parameters:\n');
fprintf('  MA Window Size: %d\n', window_size);
fprintf('  SG Order: %d\n', order);
fprintf('  SG Frame Length: %d\n', framelen);
fprintf('  WD Level: %d\n', wd_level);
fprintf('  WD Wavelet: %s\n', wd_wavelet);


% ----------------Apply noise reduction to x's-----------------
fprintf('Applying noise reduction methods...\n');
x_movmean = movmean(x_noisy, window_size); % Moving Average

% Ensure framelen is valid for sgolayfilt before applying
if mod(framelen, 2) == 0 || framelen <= order
    warning('SG Frame Length (%d) is invalid for order (%d). Using noisy data for SG part.', framelen, order);
    x_sg = x_noisy; % Fallback or handle error differently
else
    x_sg = sgolayfilt(x_noisy, order, framelen); % Savitzky-Golay
end

try
    x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet); % Wavelet Denoising
catch ME
    warning('Wavelet denoising failed with level %d, wavelet %s: %s. Using noisy data for WD part.', wd_level, wd_wavelet, ME.message);
    x_wd = x_noisy; % Fallback
end
fprintf('Noise reduction applied.\n');


% ------------------ Generate HAVOK System Data ------------------
fprintf('Generating HAVOK systems...\n');
r = TARGET_R; % Use the target r defined at the start
% Ensure getSystem is available
if exist('getSystem', 'file') ~= 2
    error('getSystem.m function not found in the MATLAB path.');
end
[V_x, A_x, B_x, xReg, U_x, E_x] = getSystem(x_original, 100, r, dt, tspan);
[V_x2, A_x2, B_x2, xReg2,U_x2, E_x2] = getSystem(x_noisy, 100, r, dt, tspan);
[V_movmean, A_x_movmean, B_x_movmean, xReg_movmean, U_movmean, E_movmean] = getSystem(x_movmean, 100, r, dt, tspan);
[V_sg, A_x_sg, B_x_sg, xReg_sg, U_sg, E_sg] = getSystem(x_sg, 100, r, dt, tspan);
[V_wd, A_x_wd, B_x_wd, xReg_wd, U_wd, E_wd] = getSystem(x_wd, 100, r, dt, tspan);
fprintf('HAVOK systems generated.\n');

% --- Check if systems are valid before proceeding ---
% Add checks here to ensure A, B, U, E matrices are not empty and have correct dimensions
% Example check:
if isempty(A_x) || isempty(A_x2) || isempty(A_x_movmean) || isempty(A_x_sg) || isempty(A_x_wd)
    error('One or more HAVOK system matrices (A) are empty. Check getSystem function or input data.');
end
% Add similar checks for B, xReg, U, E and their dimensions relative to r


% ------------------ Reconstruct and simulate all systems ------------------
fprintf('Simulating systems...\n');
L = 1:min(length(tspan), size(xReg, 1)); % Use consistent L based on original xReg size

% Define system and simulate (add error checking for each step)
try
    sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);
    [y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1)'); % Note transpose for initial condition
catch ME; warning('Simulation failed for Original: %s', '%s', ME.message); y_sim_x = nan(length(L), r-1); t_sim_x = dt*(L-1); end

try
    sys_x2 = ss(A_x2, B_x2, eye(r-1), 0*B_x2);
    [y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L, r), dt*(L-1), xReg2(1, 1:r-1)');
catch ME; warning('Simulation failed for Noisy: %s', '%s',ME.message); y_sim_x2 = nan(length(L), r-1); t_sim_x2 = dt*(L-1); end

try
    sys_x_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);
    [y_sim_x_movmean, t_sim_x_movmean] = lsim(sys_x_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1)');
catch ME; warning('Simulation failed for MovMean: %s', '%s',ME.message); y_sim_x_movmean = nan(length(L), r-1); t_sim_x_movmean = dt*(L-1); end

try
    sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);
    [y_sim_x_sg, t_sim_x_sg] = lsim(sys_x_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1)');
catch ME; warning('Simulation failed for SG: %s', '%s',ME.message); y_sim_x_sg = nan(length(L), r-1); t_sim_x_sg = dt*(L-1); end

try
    sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);
    [y_sim_x_wd, t_sim_x_wd] = lsim(sys_x_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1)');
catch ME; warning('Simulation failed for WD: %s', '%s',ME.message); y_sim_x_wd = nan(length(L), r-1); t_sim_x_wd = dt*(L-1); end

fprintf('Simulations complete.\n');


% -------------- Reconstruct x(t) for Each Method --------------
fprintf('Reconstructing x(t)...\n');
% Ensure dehankelize is available
if exist('dehankelize', 'file') ~= 2
    error('dehankelize.m function not found in the MATLAB path.');
end

% Function to safely reconstruct, handling potential errors
safe_reconstruct = @(U, E, y_sim, r_target) ...
    ( ~any(isnan(y_sim(:))) && size(U,2)>=(r_target-1) && size(E,1)>=(r_target-1) && size(E,2)>=(r_target-1) ) ...
    * dehankelize( U(:, 1:r_target-1) * E(1:r_target-1, 1:r_target-1) * y_sim.' );

x_reconstructed = safe_reconstruct(U_x, E_x, y_sim_x, r);
x_reconstructed_noisy = safe_reconstruct(U_x2, E_x2, y_sim_x2, r);
x_reconstructed_movmean = safe_reconstruct(U_movmean, E_movmean, y_sim_x_movmean, r);
x_reconstructed_sg = safe_reconstruct(U_sg, E_sg, y_sim_x_sg, r);
x_reconstructed_wd = safe_reconstruct(U_wd, E_wd, y_sim_x_wd, r);

% Check if reconstruction resulted in empty or scalar values and replace with NaNs of expected size
check_reconstruction = @(x_recon, L_len) ternary(isempty(x_recon) || isscalar(x_recon) || ~isvector(x_recon), nan(L_len, 1), x_recon(1:min(L_len, length(x_recon))));
L_len = length(L);
x_reconstructed         = check_reconstruction(x_reconstructed, L_len);
x_reconstructed_noisy   = check_reconstruction(x_reconstructed_noisy, L_len);
x_reconstructed_movmean = check_reconstruction(x_reconstructed_movmean, L_len);
x_reconstructed_sg      = check_reconstruction(x_reconstructed_sg, L_len);
x_reconstructed_wd      = check_reconstruction(x_reconstructed_wd, L_len);

fprintf('x(t) reconstructed.\n');


% -------------- Error Calculation for Reconstructed x(t) --------------
fprintf('Calculating errors...\n');
% Ensure comparison lengths are valid
compare_len = min(length(x_original), length(x_reconstructed)); % Find minimum valid length
L_comp = 1:compare_len; % Use this for comparison

% Calculate squared errors, handling potential NaNs from failed reconstructions
calc_error = @(orig, recon) ternary(any(isnan(recon(L_comp))), NaN, (orig(L_comp) - recon(L_comp)).^2);

error_x           = calc_error(x_original, x_reconstructed);
error_x_noisy     = calc_error(x_original, x_reconstructed_noisy);
error_x_movmean   = calc_error(x_original, x_reconstructed_movmean);
error_x_sg        = calc_error(x_original, x_reconstructed_sg);
error_x_wd        = calc_error(x_original, x_reconstructed_wd);

fprintf('Errors calculated.\n');


% -------------- Plot Results for Comparison --------------
fprintf('Generating plots...\n');

% Plot the noisy delay emebdded attractors
figure('Name', sprintf('Noisy Attractor (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
subplot(2, 2, 1);
plot3(V_x(L_comp,1), V_x(L_comp,2), V_x(L_comp,3), 'b', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Original Attractor');
grid on; view(3);
subplot(2, 2, 2);
plot3(V_x2(L_comp,1), V_x2(L_comp,2), V_x2(L_comp,3), 'r', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Noisy Attractor');
grid on; view(3);
subplot(2, 2, 3);
plot3(V_sg(L_comp,1), V_sg(L_comp,2), V_sg(L_comp,3), 'm', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Savitzky-Golay Attractor');
grid on; view(3);
subplot(2, 2, 4);
plot3(V_wd(L_comp,1), V_wd(L_comp,2), V_wd(L_comp,3), 'k', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Wavelet Attractor');
grid on; view(3);

%Plot the reconstructed embedded attractors
figure('Name', sprintf('Reconstructed Attractor (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
subplot(2, 2, 1);
plot3(y_sim_x(L_comp,1), y_sim_x(L_comp,2), y_sim_x(L_comp,3), 'b', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Original Attractor');
grid on; view(3);
subplot(2, 2, 2);
plot3(y_sim_x2(L_comp,1), y_sim_x2(L_comp,2), y_sim_x2(L_comp,3), 'r', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Noisy Attractor');
grid on; view(3);
subplot(2, 2, 3);
plot3(y_sim_x_sg(L_comp,1), y_sim_x_sg(L_comp,2), y_sim_x_sg(L_comp,3), 'm', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Savitzky-Golay Attractor');
grid on; view(3);
subplot(2, 2, 4);
plot3(y_sim_x_wd(L_comp,1), y_sim_x_wd(L_comp,2), y_sim_x_wd(L_comp,3), 'k', 'LineWidth', 1.2);
xlabel('v1'); ylabel('v2'); zlabel('v3');
title('Wavelet Attractor');
grid on; view(3);



figure('Name', sprintf('Reconstruction Comparison (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
hold on;
plot(t_sim_x, x_reconstructed(L_comp), 'b', 'LineWidth', 1.2);
plot(t_sim_x2, x_reconstructed_noisy(L_comp), 'r', 'LineWidth', 1.2);
plot(t_sim_x_movmean, x_reconstructed_movmean(L_comp), 'g', 'LineWidth', 1.2);
plot(t_sim_x_sg, x_reconstructed_sg(L_comp), 'm', 'LineWidth', 1.2);
plot(t_sim_x_wd, x_reconstructed_wd(L_comp), 'c', 'LineWidth', 1.2); % Cyan for Wavelet
legend('HAVOK Original', 'Noisy HAVOK', 'Moving Avg', 'Savitzky-Golay', 'Wavelet');
xlabel('Time'); ylabel('x(t)');
title(sprintf('Reconstruction of x(t) (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
grid on;
hold off;

% ... (Keep your other plotting sections for Attractors, Error Comparison, Bar Charts etc.) ...
% Make sure to use L_comp for indexing where appropriate, especially for error plots
% and comparisons involving x_original.

% Example modification for Error Plot:
figure('Name', sprintf('Error Comparison (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
hold on;
plot(tspan(L_comp), error_x, 'b', 'LineWidth', 1);
plot(tspan(L_comp), error_x_noisy, 'r', 'LineWidth', 1);
plot(tspan(L_comp), error_x_movmean, 'g', 'LineWidth', 1);
plot(tspan(L_comp), error_x_sg, 'm', 'LineWidth', 1);
plot(tspan(L_comp), error_x_wd, 'c', 'LineWidth', 1); % Cyan for Wavelet
title('Squared Error Comparison');
xlabel('Time'); ylabel('Squared Error');
legend('Original','Noisy' , 'Moving Average', 'Savitzky-Golay', 'Wavelet');
grid on;
hold off;


% Example modification for Mean Error Bar Chart:
figure('Name', sprintf('Mean Error Bar Chart (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
errorMeans = [mean(error_x, 'omitnan'), mean(error_x_noisy, 'omitnan'), mean(error_x_movmean, 'omitnan'), ...
             mean(error_x_sg, 'omitnan'), mean(error_x_wd, 'omitnan')];
methods = {'HAVOK Original', 'Noisy HAVOK', 'Moving Average', 'Savitzky-Golay', 'Wavelet'};
bar(errorMeans);
set(gca, 'XTickLabel', methods, 'XTick',1:numel(methods), 'XTickLabelRotation', 15);
ylabel('Mean of Squared Errors');
title(sprintf('Mean Errors for Each Method (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
grid on;

% --- Add the rest of your plotting code here ---
% Remember to use L_comp or appropriate indices when comparing with x_original
% or plotting errors. For attractor plots using simulated y_sim variables,
% you can still use L.

%Plot separate figures for each noise reduction method

% Moving Average

Rs = [1, int8(r/2), r-1];
figure('Name', sprintf('Reconstructed Vs with MA (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_movmean(L,i), 'g', 'LineWidth', 1);
    title(['Moving Average: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Moving Average');
end

% Savitzky-Golay
figure('Name', sprintf('Reconstructed Vs with SG (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_sg(L,i), 'm', 'LineWidth', 1);
    title(['Savitzky-Golay: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Savitzky-Golay');
end

% Wavelet Denoising
figure('Name', sprintf('Reconstructed Vs with WD (r=%d, var=%g)', TARGET_R, TARGET_VARIANCE));
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_wd(L,i), 'k', 'LineWidth', 1);
    title(['Wavelet: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Wavelet');
end


fprintf('Analysis and plotting complete.\n');
% ------------------ End of Script ------------------

% Helper for conditional assignment (ternary operator)
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end