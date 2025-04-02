% run_all_optimizations.m
clear; clc; close all;

fprintf('Starting precomputation of optimal denoising parameters...\n');

% Load the original data
try
    load('Data/lorenzData.mat', 'sol', 't', 'dt'); % Contains 'sol', 't', 'dt'
    x_original = sol(:,1);
    tspan = dt:dt:t(end); % Use full available timespan from loaded data
    fprintf('Loaded Lorenz data successfully.\n');
catch ME
    error('Failed to load Data/lorenzData.mat. Ensure the file exists and is in the correct path. Error: %s', ME.message);
end

% Define parameter space for optimization
r_values = [4, 5, 6, 7, 8, 9, 10];
variance_values = [1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001];

% Define parameter ranges to be tested by the optimization function
% These can be adjusted as needed
paramRanges.window_sizes = [3, 5, 7, 9, 11, 13, 15];
paramRanges.sg_orders    = [2, 3, 4];
paramRanges.sg_framelens = [11:4:111]; % Must be odd and > order
paramRanges.wd_levels    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
paramRanges.wd_wavelets  = {'haar', 'db4', 'sym4', 'coif3', 'bior3.5', 'dmey'};

% Initialize structure to store results
% Using nested structs: optimization_results.r{r_val}.var{var_key_str}
optimization_results = struct();

% Seed random number generator for reproducibility of noise
rng(1);

total_combinations = length(r_values) * length(variance_values);
current_combination = 0;

% Loop through all combinations
for r_val = r_values
    for var_val = variance_values
        current_combination = current_combination + 1;
        fprintf('\n============================================================\n');
        fprintf('Processing Combination %d/%d: r = %d, variance = %g\n', ...
                current_combination, total_combinations, r_val, var_val);
        fprintf('============================================================\n');

        % Generate noisy data for this variance
        x_noisy = x_original + sqrt(var_val) * randn(size(x_original)); % Use sqrt(variance) for std dev

        % Call the optimization function
        % Ensure findOptimalDenoiseParams is on the MATLAB path
        try
            bestParams = findOptimalDenoiseParams(x_original, x_noisy, r_val, dt, tspan, paramRanges);

            % Create a valid field name for the variance (replace '.' with 'p')
            var_key_str = sprintf('var%g', var_val);
            var_key_str = strrep(var_key_str, '.', 'p');
            var_key_str = strrep(var_key_str, '-', 'neg'); % Handle negative if ever needed

            % Store the results (only parameters)
            optimization_results.(sprintf('r%d', r_val)).(var_key_str) = bestParams;

            fprintf('Stored results for r=%d, variance=%g\n', r_val, var_val);

        catch ME
             fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
             fprintf('ERROR processing combination r=%d, variance=%g\n', r_val, var_val);
             fprintf('Error message: %s\n', ME.message);
             fprintf('Occurred in %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
             fprintf('Skipping this combination.\n');
             fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        end
    end
end

% Save the results structure to a .mat file
results_filename = 'optimal_denoise_params.mat';
try
    save(results_filename, 'optimization_results');
    fprintf('\n\nPrecomputation complete. Optimal parameters saved to %s\n', results_filename);
catch ME
    error('Failed to save results to %s. Error: %s', results_filename, ME.message);
end

% --- Make sure helper functions are available ---
% Ensure getSystem.m and dehankelize.m are in the MATLAB path or current directory.
% You might need to add their directory:
% addpath('path/to/helper/functions');