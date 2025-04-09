% compare_denoising_across_variances.m
clear; clc; close all;

% -------- USER DEFINED PARAMETERS --------
TARGET_R = 7; % Set the desired HAVOK rank for COMPARISON ACROSS VARIANCES
variance_values = [25, 10, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.00001, 0.0000001]; % Variances to test
% variance_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]; % Example: Lower variances

results_filename = 'Data/optimal_denoise_params.mat'; % File with precomputed params
run_optimization_if_missing = true; % Set to true to run optimizer if params not found
save_newly_computed_params = true; % Set to true to save parameters if computed on-the-fly
% -----------------------------------------

fprintf('Starting comparison for r = %d across %d variances.\n', TARGET_R, length(variance_values));

% --- Load the original data ---
try
    load('Data/lorenzData.mat', 'sol', 't', 'dt'); % Contains 'sol', 't', 'dt'
    x_original = sol(:,1);
    tspan = dt:dt:t(end);
    fprintf('Loaded Lorenz data.\n');
catch ME
    error('Failed to load Data/lorenzData.mat: %s', ME.message);
end

% --- Ensure Helper Functions Exist ---
if exist('getSystem', 'file') ~= 2
    error('getSystem.m function not found in the MATLAB path.');
end
if exist('dehankelize', 'file') ~= 2
    error('dehankelize.m function not found in the MATLAB path.');
end
if exist('findOptimalDenoiseParams', 'file') ~= 2 && run_optimization_if_missing
    warning('findOptimalDenoiseParams.m function not found. On-the-fly optimization will fail if needed.');
end

% --- Pre-calculate Original System and Reconstruction (Noise-Free Baseline) ---
fprintf('Generating noise-free baseline HAVOK system (r=%d)...\n', TARGET_R);
try
    [V_x, A_x, B_x, xReg, U_x, E_x] = getSystem(x_original, 100, TARGET_R, dt, tspan);
    sys_x = ss(A_x, B_x, eye(TARGET_R-1), 0*B_x);
    sim_len = min(length(tspan), size(xReg, 1)); % Determine length for simulation/comparison
    L_comp = 1:sim_len;
    [y_sim_x, ~] = lsim(sys_x, xReg(L_comp, TARGET_R), dt*(L_comp-1), xReg(1, 1:TARGET_R-1)');
    x_reconstructed_original = dehankelize( U_x(:, 1:TARGET_R-1) * E_x(1:TARGET_R-1, 1:TARGET_R-1) * y_sim_x.' );
    x_reconstructed_original = x_reconstructed_original(1:min(length(L_comp), length(x_reconstructed_original))); % Ensure length consistency
    if length(x_reconstructed_original) < length(L_comp)
       warning('Original reconstruction shorter than expected simulation length.');
       % Pad with NaN or handle as appropriate if needed, L_comp might need adjustment
       x_reconstructed_original = [x_reconstructed_original; nan(length(L_comp)-length(x_reconstructed_original),1)];
    end
    mse_original = mean((x_original(L_comp) - x_reconstructed_original(L_comp)).^2, 'omitnan');
    fprintf('Noise-free baseline MSE: %g\n', mse_original);
catch ME
    error('Failed to generate or simulate the noise-free baseline system: %s', ME.message);
end


% --- Initialize Storage for MSE Results ---
mse_noisy   = nan(1, length(variance_values));
mse_movmean = nan(1, length(variance_values));
mse_sg      = nan(1, length(variance_values));
mse_wd      = nan(1, length(variance_values));

% --- Load Optimization Results File Once ---
optimization_results = struct(); % Initialize empty
if exist(results_filename, 'file')
    fprintf('Loading precomputed parameters from %s...\n', results_filename);
    load(results_filename, 'optimization_results'); % Loads 'optimization_results' struct
    fprintf('Precomputed parameters loaded.\n');
else
    fprintf('Precomputed parameter file %s not found.\n', results_filename);
    if ~run_optimization_if_missing
       error('Parameter file not found and on-the-fly optimization is disabled. Cannot proceed.');
    end
end

% --- Loop Through Each Variance ---
needs_saving = false; % Flag to track if new results were computed and need saving

for idx = 1:length(variance_values)
    current_variance = variance_values(idx);
    fprintf('\n--- Processing Variance = %g (%d/%d) ---\n', current_variance, idx, length(variance_values));

    % --- Add Gaussian Noise ---
    rng('default'); % Reset RNG for consistent noise generation for THIS variance level across runs
    rng(idx);       % Use index for seed to get different noise patterns for different variances
    x_noisy = x_original + sqrt(current_variance) * randn(size(x_original));
    fprintf('Added Gaussian noise with variance %g.\n', current_variance);

    % --- Get Optimal Denoising Parameters for Current R and Variance ---
    params_found = false;
    bestParams = struct();
    r_field = sprintf('r%d', TARGET_R);
    var_field = sprintf('var%g', current_variance);
    var_field = strrep(var_field, '.', 'p'); % Replace '.' with 'p' for valid field name
    var_field = strrep(var_field, '-', 'neg'); % Handle negative exponents if they were used

    if isfield(optimization_results, r_field) && isfield(optimization_results.(r_field), var_field)
        bestParams = optimization_results.(r_field).(var_field);
        if isstruct(bestParams) && ~isempty(fieldnames(bestParams)) % Basic check
             params_found = true;
             fprintf('Found precomputed parameters.\n');
        else
            fprintf('Precomputed parameter entry found but invalid for r=%d, var=%g.\n', TARGET_R, current_variance);
        end
    else
        fprintf('Precomputed parameters NOT found for r=%d, var=%g.\n', TARGET_R, current_variance);
    end

    % --- Run Optimization if Parameters Not Found ---
    if ~params_found
        if run_optimization_if_missing && exist('findOptimalDenoiseParams', 'file') == 2
            fprintf('Running optimization on-the-fly...\n');

            % Define parameter ranges FOR ON-THE-FLY optimization
            paramRanges.window_sizes = [3, 5, 7, 9, 11, 13, 15];
            paramRanges.sg_orders    = [2, 3, 4];
            paramRanges.sg_framelens = [11:4:111];
            paramRanges.wd_levels    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
            paramRanges.wd_wavelets  = {'haar', 'db4', 'sym4', 'coif3', 'bior3.5', 'dmey'};

            try
                bestParams = findOptimalDenoiseParams(x_original, x_noisy, TARGET_R, dt, tspan, paramRanges);
                params_found = true;
                fprintf('On-the-fly optimization complete.\n');

                % Store newly computed parameters
                 if ~isfield(optimization_results, r_field)
                    optimization_results.(r_field) = struct();
                 end
                 optimization_results.(r_field).(var_field) = bestParams;
                 needs_saving = true; % Mark that we have new results to potentially save

            catch ME
                warning('Error during on-the-fly optimization for var=%g: %s\nSkipping this variance.', current_variance, ME.message);
                continue; % Skip to the next variance value
            end
        else
            warning('Required parameters not found and cannot run optimization. Skipping variance %g.', current_variance);
            continue; % Skip to the next variance value
        end
    end

    % --- Assign Parameters ---
    if ~params_found || isempty(fieldnames(bestParams))
        warning('Failed to obtain valid denoising parameters for variance %g. Skipping.', current_variance);
        continue; % Skip to the next variance value
    end

    window_size = bestParams.ma_window_size;
    order       = bestParams.sg_order;
    framelen    = bestParams.sg_framelen;
    wd_level    = bestParams.wd_level;
    wd_wavelet  = bestParams.wd_wavelet;

    fprintf('Using Parameters: MA Win=%d, SG Ord=%d, SG Frame=%d, WD Lvl=%d, WD Wav=%s\n', ...
            window_size, order, framelen, wd_level, wd_wavelet);

    % --- Apply Noise Reduction ---
    fprintf('Applying noise reduction methods...\n');
    x_movmean = movmean(x_noisy, window_size);

    if mod(framelen, 2) == 0 || framelen <= order
        warning('SG Frame Length (%d) is invalid for order (%d) at var=%g. Using noisy data for SG part.', framelen, order, current_variance);
        x_sg = x_noisy;
    else
        x_sg = sgolayfilt(x_noisy, order, framelen);
    end

    try
        x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet);
    catch ME
        warning('Wavelet denoising failed for var=%g (level %d, wavelet %s): %s. Using noisy data for WD part.', current_variance, wd_level, wd_wavelet, ME.message);
        x_wd = x_noisy;
    end
    fprintf('Noise reduction applied.\n');

    % --- Generate HAVOK Systems for Noisy and Denoised Data ---
    fprintf('Generating HAVOK systems...\n');
    try
        [~, A_x2, B_x2, xReg2, U_x2, E_x2] = getSystem(x_noisy, 100, TARGET_R, dt, tspan);
        [~, A_movmean, B_movmean, xReg_movmean, U_movmean, E_movmean] = getSystem(x_movmean, 100, TARGET_R, dt, tspan);
        [~, A_sg, B_sg, xReg_sg, U_sg, E_sg] = getSystem(x_sg, 100, TARGET_R, dt, tspan);
        [~, A_wd, B_wd, xReg_wd, U_wd, E_wd] = getSystem(x_wd, 100, TARGET_R, dt, tspan);
        fprintf('HAVOK systems generated.\n');
    catch ME
         warning('Failed to generate one or more HAVOK systems for var=%g: %s Skipping.', current_variance, ME.message);
         continue; % Skip to next variance
    end


    % --- Simulate Systems & Reconstruct x(t) ---
    fprintf('Simulating systems and reconstructing x(t)...\n');
    % Define a safe reconstruction function inline
    safe_reconstruct = @(U, E, y_sim, r_target) ...
        ternary( (~any(isnan(y_sim(:))) && size(U,2)>=(r_target-1) && size(E,1)>=(r_target-1) && size(E,2)>=(r_target-1) && r_target>1), ...
                 dehankelize( U(:, 1:r_target-1) * E(1:r_target-1, 1:r_target-1) * y_sim.' ), ...
                 nan(size(x_original,1), 1) ); % Return NaNs on failure

    % Simulate and Reconstruct for each method
    try
        sys_x2 = ss(A_x2, B_x2, eye(TARGET_R-1), 0*B_x2);
        [y_sim_x2, ~] = lsim(sys_x2, xReg2(L_comp, TARGET_R), dt*(L_comp-1), xReg2(1, 1:TARGET_R-1)');
        x_reconstructed_noisy = safe_reconstruct(U_x2, E_x2, y_sim_x2, TARGET_R);
    catch ME; warning('Sim/Recon failed for Noisy (var=%g): %s', current_variance, ME.message); x_reconstructed_noisy = nan(length(L_comp), 1); end

    try
        sys_movmean = ss(A_movmean, B_movmean, eye(TARGET_R-1), 0*B_movmean);
        [y_sim_movmean, ~] = lsim(sys_movmean, xReg_movmean(L_comp, TARGET_R), dt*(L_comp-1), xReg_movmean(1, 1:TARGET_R-1)');
        x_reconstructed_movmean = safe_reconstruct(U_movmean, E_movmean, y_sim_movmean, TARGET_R);
    catch ME; warning('Sim/Recon failed for MovMean (var=%g): %s', current_variance, ME.message); x_reconstructed_movmean = nan(length(L_comp), 1); end

    try
        sys_sg = ss(A_sg, B_sg, eye(TARGET_R-1), 0*B_sg);
        [y_sim_sg, ~] = lsim(sys_sg, xReg_sg(L_comp, TARGET_R), dt*(L_comp-1), xReg_sg(1, 1:TARGET_R-1)');
        x_reconstructed_sg = safe_reconstruct(U_sg, E_sg, y_sim_sg, TARGET_R);
    catch ME; warning('Sim/Recon failed for SG (var=%g): %s', current_variance, ME.message); x_reconstructed_sg = nan(length(L_comp), 1); end

    try
        sys_wd = ss(A_wd, B_wd, eye(TARGET_R-1), 0*B_wd);
        [y_sim_wd, ~] = lsim(sys_wd, xReg_wd(L_comp, TARGET_R), dt*(L_comp-1), xReg_wd(1, 1:TARGET_R-1)');
        x_reconstructed_wd = safe_reconstruct(U_wd, E_wd, y_sim_wd, TARGET_R);
    catch ME; warning('Sim/Recon failed for WD (var=%g): %s', current_variance, ME.message); x_reconstructed_wd = nan(length(L_comp), 1); end

    fprintf('Simulations and reconstructions complete.\n');

    % --- Calculate and Store Mean Squared Error ---
    fprintf('Calculating and storing MSE...\n');
    calc_mse = @(orig, recon, len) mean((orig(1:len) - recon(1:min(len, length(recon)))).^2, 'omitnan');

    current_len = length(L_comp); % Use the length determined from the original simulation

    mse_noisy(idx)   = calc_mse(x_original, x_reconstructed_noisy, current_len);
    mse_movmean(idx) = calc_mse(x_original, x_reconstructed_movmean, current_len);
    mse_sg(idx)      = calc_mse(x_original, x_reconstructed_sg, current_len);
    mse_wd(idx)      = calc_mse(x_original, x_reconstructed_wd, current_len);

    fprintf('MSEs: Noisy=%.4g, MA=%.4g, SG=%.4g, WD=%.4g\n', mse_noisy(idx), mse_movmean(idx), mse_sg(idx), mse_wd(idx));

end % End of variance loop

% --- Optional: Save Updated Parameter File ---
if needs_saving && save_newly_computed_params
    save_confirm = input(sprintf('Save newly computed parameters to %s? (y/n): ', results_filename), 's');
    if lower(save_confirm) == 'y'
        try
            save(results_filename, 'optimization_results');
            fprintf('Saved updated parameters to %s\n', results_filename);
        catch ME
            warning(ME.identifier, 'Failed to save updated parameters: %s', ME.message);
        end
    else
        fprintf('Skipping save of newly computed parameters.\n');
    end
end


% --- Plot Final Comparison Figure ---
fprintf('\nGenerating final comparison plot...\n');
figure('Name', sprintf('Denoising Method Comparison (r=%d)', TARGET_R), 'Position', [100, 100, 800, 600]);
hold on;

% Plot MSE for each method against variance
plot(variance_values, mse_noisy, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(variance_values, mse_movmean, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(variance_values, mse_sg, 'm-^', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(variance_values, mse_wd, 'c-d', 'LineWidth', 1.5, 'MarkerSize', 6);

% Plot baseline original MSE (should be near constant)
plot(variance_values, repmat(mse_original, size(variance_values)), 'b--', 'LineWidth', 1.5);

% --- Plot Customization ---
set(gca, 'XScale', 'log'); % Use log scale for variance x-axis
set(gca, 'YScale', 'log'); % Use log scale for MSE y-axis (often useful for errors)
% set(gca, 'YScale', 'linear'); % Use linear scale if preferred

xlabel('Noise Variance (\sigma^2)');
ylabel('Mean Squared Error (MSE) of Reconstruction');
title(sprintf('HAVOK Reconstruction Error vs. Noise Variance (r = %d)', TARGET_R));
legend('Noisy Reconstruction', 'Moving Average Denoised', 'Savitzky-Golay Denoised', 'Wavelet Denoised', 'Original (Noise-Free)', 'Location', 'best');
grid on;
set(gca, 'XDir', 'reverse'); % Show highest variance on the left, typical for noise plots
hold off;

fprintf('Analysis and plotting complete.\n');
% ------------------ End of Script ------------------

% Helper for conditional assignment (ternary operator) - Keep this function available
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end