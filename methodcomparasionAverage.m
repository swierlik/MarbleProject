% COMPARENOISEMETHODSMULTIPLEREALIZATIONS Compares noise reduction methods across
% multiple noise realizations at a specified noise level
%
% Inputs:
%   noise_level - The standard deviation of noise to add
%   num_realizations - Number of different noise realizations to test (default: 25)
%
% Example:
%   compareNoiseMethodsMultipleRealizations(0.1, 25);

% Load the original deterministic Lorenz data
load('Data/lorenzData.mat');

noise_level = 0.0001; % Standard deviation of noise to add
num_realizations = 5; % Number of different noise realizations to test

% Define parameters for noise reduction methods
window_size = 3; % Moving Average
sg_order = 2; % Savitzky-Golay Order
sg_framelen = 45; % Savitzky-Golay Frame Length
wd_level = 1; % Wavelet Denoising Level
wd_wavelet = 'bior3.5'; % Wavelet Denoising Wavelet

% Initialize arrays to store MSE for each method and realization
mse_noisy = zeros(num_realizations, 1);
mse_movmean = zeros(num_realizations, 1);
mse_sg = zeros(num_realizations, 1);
mse_wd = zeros(num_realizations, 1);

% Define a maximum cap for MSE values to prevent extreme outliers
mse_cap = 1000; % Adjust this value based on your typical MSE range

% Extract original x data
x_original = sol(:,1);

% Loop through each noise realization
for i = 1:num_realizations
    fprintf('Processing noise realization %d of %d\n', i, num_realizations);
    
    % Set random seed based on realization number for reproducibility
    % but different noise pattern for each realization
    rng(i); 
    
    % Generate noisy data with specified noise level
    sol_noisy = sol;
    sol_noisy(:,1) = sol(:,1) + noise_level * randn(size(sol,1), 1);
    
    % Extract noisy x data
    x_noisy = sol_noisy(:,1);
    
    % Apply noise reduction methods
    x_movmean = movmean(x_noisy, window_size);
    x_sg = sgolayfilt(x_noisy, sg_order, sg_framelen);
    x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet);
    
    % Generate HAVOK systems for each method
    r = 15; % Rank for HAVOK model
    
    % Process each dataset - wrap in try-catch to handle potential numerical issues
    try
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
        
        % Simulate systems (wrap in try-catch for each simulation)
        try
            [y_sim_orig, t_sim_orig] = lsim(sys_orig, xReg_orig(L, r), dt*(L-1), xReg_orig(1, 1:r-1));
            [y_sim_noisy, t_sim_noisy] = lsim(sys_noisy, xReg_noisy(L, r), dt*(L-1), xReg_noisy(1, 1:r-1));
            
            % Ensure correct dimensions for reconstruction
            U_orig = U_orig(:, 1:r-1);
            U_noisy = U_noisy(:, 1:r-1);
            
            E_orig = E_orig(1:r-1, 1:r-1);
            E_noisy = E_noisy(1:r-1, 1:r-1);
            
            % Reconstruct signals
            x_reconstructed_orig = U_orig * E_orig * y_sim_orig.';
            x_reconstructed_noisy = U_noisy * E_noisy * y_sim_noisy.';
            
            % Calculate MSE for noisy case
            ref_x = x_original(L);
            min_length = min(length(ref_x), size(x_reconstructed_noisy, 2));
            ref_x_trimmed = ref_x(1:min_length)';
            x_reconstructed_noisy_trimmed = x_reconstructed_noisy(1, 1:min_length);
            mse_noisy(i) = mean((ref_x_trimmed - x_reconstructed_noisy_trimmed).^2);
        catch
            fprintf('Error simulating noisy system in realization %d. Using capped value.\n', i);
            mse_noisy(i) = mse_cap;
        end
        
        % Moving Average - separate try-catch for each method
        try
            [y_sim_movmean, t_sim_movmean] = lsim(sys_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));
            U_movmean = U_movmean(:, 1:r-1);
            E_movmean = E_movmean(1:r-1, 1:r-1);
            x_reconstructed_movmean = U_movmean * E_movmean * y_sim_movmean.';
            
            ref_x = x_original(L);
            min_length = min(length(ref_x), size(x_reconstructed_movmean, 2));
            ref_x_trimmed = ref_x(1:min_length)';
            x_reconstructed_movmean_trimmed = x_reconstructed_movmean(1, 1:min_length);
            
            mse_value = mean((ref_x_trimmed - x_reconstructed_movmean_trimmed).^2);
            % Check for problematic values and cap them
            if isnan(mse_value) || isinf(mse_value) || mse_value > mse_cap
                mse_movmean(i) = mse_cap;
                fprintf('Capping problematic moving average MSE in realization %d\n', i);
            else
                mse_movmean(i) = mse_value;
            end
        catch
            fprintf('Error simulating moving average in realization %d. Using capped value.\n', i);
            mse_movmean(i) = mse_cap;
        end
        
        % Savitzky-Golay
        try
            [y_sim_sg, t_sim_sg] = lsim(sys_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));
            U_sg = U_sg(:, 1:r-1);
            E_sg = E_sg(1:r-1, 1:r-1);
            x_reconstructed_sg = U_sg * E_sg * y_sim_sg.';
            
            ref_x = x_original(L);
            min_length = min(length(ref_x), size(x_reconstructed_sg, 2));
            ref_x_trimmed = ref_x(1:min_length)';
            x_reconstructed_sg_trimmed = x_reconstructed_sg(1, 1:min_length);
            
            mse_value = mean((ref_x_trimmed - x_reconstructed_sg_trimmed).^2);
            if isnan(mse_value) || isinf(mse_value) || mse_value > mse_cap
                mse_sg(i) = mse_cap;
                fprintf('Capping problematic Savitzky-Golay MSE in realization %d\n', i);
            else
                mse_sg(i) = mse_value;
            end
        catch
            fprintf('Error simulating Savitzky-Golay in realization %d. Using capped value.\n', i);
            mse_sg(i) = mse_cap;
        end
        
        % Wavelet
        try
            [y_sim_wd, t_sim_wd] = lsim(sys_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));
            U_wd = U_wd(:, 1:r-1);
            E_wd = E_wd(1:r-1, 1:r-1);
            x_reconstructed_wd = U_wd * E_wd * y_sim_wd.';
            
            ref_x = x_original(L);
            min_length = min(length(ref_x), size(x_reconstructed_wd, 2));
            ref_x_trimmed = ref_x(1:min_length)';
            x_reconstructed_wd_trimmed = x_reconstructed_wd(1, 1:min_length);
            
            mse_value = mean((ref_x_trimmed - x_reconstructed_wd_trimmed).^2);
            if isnan(mse_value) || isinf(mse_value) || mse_value > mse_cap
                mse_wd(i) = mse_cap;
                fprintf('Capping problematic Wavelet MSE in realization %d\n', i);
            else
                mse_wd(i) = mse_value;
            end
        catch
            fprintf('Error simulating Wavelet in realization %d. Using capped value.\n', i);
            mse_wd(i) = mse_cap;
        end
    catch
        fprintf('Error in processing realization %d. Setting all methods to capped value.\n', i);
        mse_noisy(i) = mse_cap;
        mse_movmean(i) = mse_cap;
        mse_sg(i) = mse_cap;
        mse_wd(i) = mse_cap;
    end
end

% Calculate statistics (with handling for potential remaining Inf/NaN values)
mse_noisy(isinf(mse_noisy) | isnan(mse_noisy)) = mse_cap;
mse_movmean(isinf(mse_movmean) | isnan(mse_movmean)) = mse_cap;
mse_sg(isinf(mse_sg) | isnan(mse_sg)) = mse_cap;
mse_wd(isinf(mse_wd) | isnan(mse_wd)) = mse_cap;

mean_mse = [mean(mse_noisy), mean(mse_movmean), mean(mse_sg), mean(mse_wd)];
std_mse = [std(mse_noisy), std(mse_movmean), std(mse_sg), std(mse_wd)];

% Check if any method is at cap for all realizations and report
methods = {'No Filtering', 'Moving Average', 'Savitzky-Golay', 'Wavelet'};
all_mse = {mse_noisy, mse_movmean, mse_sg, mse_wd};
for i = 1:4
    if all(all_mse{i} >= mse_cap * 0.99)
        fprintf('Warning: %s method reached cap value for all realizations\n', methods{i});
    end
end

% Create bar plot
figure;
bar_h = bar(mean_mse);
hold on;

% Add error bars
errorbar(1:4, mean_mse, std_mse, 'k.', 'LineWidth', 1.5);

% Customize plot
set(gca, 'XTickLabel', methods);
xlabel('Noise Reduction Method', 'FontSize', 14);
ylabel('Mean Square Error', 'FontSize', 14);
title(['MSE Comparison Across ', num2str(num_realizations), ' Noise Realizations (η = ', num2str(noise_level), ')'], 'FontSize', 16);
grid on;
ylim([0, mse_cap]);

% Add text showing improvement percentages relative to no filtering
for i = 2:4
    improvement = (1 - mean_mse(i)/mean_mse(1)) * 100;
    % Only show improvement text if it's not at the cap
    if mean_mse(i) < mse_cap * 0.99
        text(i, mean_mse(i)/2, [num2str(improvement, '%.1f'), '% better'], ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontWeight', 'bold');
    else
        text(i, mean_mse(i)/2, 'Unstable', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontWeight', 'bold');
    end
end

% Optional: Log-scale for better visualization if the differences are large
if max(mean_mse)/min(mean_mse) > 10
    set(gca, 'YScale', 'log');
    ylabel('Mean Square Error (log scale)', 'FontSize', 14);
end

% Create scatter plot to see distribution of results
figure;
colors = ['r', 'g', 'b', 'k'];
for i = 1:4
    switch i
        case 1
            data = mse_noisy;
        case 2
            data = mse_movmean;
        case 3
            data = mse_sg;
        case 4
            data = mse_wd;
    end
    
    % Create jittered x positions
    x = i + (rand(size(data)) - 0.5) * 0.3;
    
    % Plot points
    scatter(x, data, 50, colors(i), 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    
    % Add mean and std lines
    line([i-0.2, i+0.2], [mean(data), mean(data)], 'Color', colors(i), 'LineWidth', 2);
    line([i, i], [mean(data)-std(data), mean(data)+std(data)], 'Color', colors(i), 'LineWidth', 2);
end

% Customize scatter plot
set(gca, 'XTick', 1:4, 'XTickLabel', methods);
xlabel('Noise Reduction Method', 'FontSize', 14);
ylabel('Mean Square Error', 'FontSize', 14);
title(['MSE Distribution Across ', num2str(num_realizations), ' Noise Realizations (η = ', num2str(noise_level), ')'], 'FontSize', 16);
grid on;

% Optional: Log-scale for better visualization if the differences are large
if max([max(mse_noisy), max(mse_movmean), max(mse_sg), max(mse_wd)]) / ...
   min([min(mse_noisy), min(mse_movmean), min(mse_sg), min(mse_wd)]) > 10
    set(gca, 'YScale', 'log');
    ylabel('Mean Square Error (log scale)', 'FontSize', 14);
end

% Print statistics to console
fprintf('\n--- Results Summary (Noise Level: %.6g) ---\n', noise_level);
fprintf('Method          | Mean MSE      | Std Dev       | Improvement\n');
fprintf('----------------|---------------|---------------|-------------\n');
fprintf('No Filtering    | %.6e | %.6e | Baseline\n', mean_mse(1), std_mse(1));

for i = 2:4
    improvement = (1 - mean_mse(i)/mean_mse(1)) * 100;
    if mean_mse(i) >= mse_cap * 0.99
        fprintf('%s | %.6e | %.6e | Unstable (capped)\n', ...
            pad(methods{i}, 15), mean_mse(i), std_mse(i));
    else
        fprintf('%s | %.6e | %.6e | %.2f%%\n', ...
            pad(methods{i}, 15), mean_mse(i), std_mse(i), improvement);
    end
end

% Save results
results = struct();
results.noise_level = noise_level;
results.num_realizations = num_realizations;
results.mse_noisy = mse_noisy;
results.mse_movmean = mse_movmean;
results.mse_sg = mse_sg;
results.mse_wd = mse_wd;
results.mean_mse = mean_mse;
results.std_mse = std_mse;
results.mse_cap = mse_cap;

save(['noise_comparison_eta_', num2str(noise_level), '_capped.mat'], 'results');
fprintf('\nResults saved to noise_comparison_eta_%.6g_capped.mat\n', noise_level);