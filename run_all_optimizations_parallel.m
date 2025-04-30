% run_all_optimizations_parallel.m
% Runs optimizations for multiple parameter combinations in parallel
clear; clc; close all;

fprintf('Starting precomputation of optimal denoising parameters with parallel processing...\n');

clock_start = tic; % Start timer
% % Load the original data
% try
%     load('Data/lorenzData.mat', 'sol', 't', 'dt'); % Contains 'sol', 't', 'dt'
%     x_original = sol(:,1);
%     tspan = dt:dt:t(end); % Use full available timespan from loaded data
%     fprintf('Loaded Lorenz data successfully.\n');
% catch ME
%     error('Failed to load Data/lorenzData.mat. Ensure the file exists and is in the correct path. Error: %s', ME.message);
% end

% Define parameter space for optimization
r_values = [4, 5, 6, 7, 8, 9, 10];
variance_values = [100, 10, 1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.00001, 0.0000001, 0.00000001];

% Number of runs for averaging in each parameter combination
num_runs = 10; % Adjust based on computational resources

%Lorenz data initial conditions
initial_conditions = [-8, 8, 27];
dt=0.001;
tspan = [dt:dt:50]; % Time span for the Lorenz system

% Define parameter ranges to be tested by the optimization function
paramRanges = struct();
paramRanges.window_sizes = [3, 5, 7, 9, 11, 13, 15];
paramRanges.sg_orders    = [2, 3, 4];
paramRanges.sg_framelens = [11:4:111]; % Must be odd and > order
paramRanges.wd_levels    = 1:2:20;
paramRanges.wd_wavelets  = {'haar', 'db4', 'sym4', 'coif3', 'bior3.5', 'dmey'};

% Initialize parallel pool if not already running
if isempty(gcp('nocreate'))
    try
        parpool('local'); % Start with default number of workers
        fprintf('Parallel pool started successfully.\n');
    catch ME
        warning(ME.identifier, 'Could not start parallel pool: %s\n', ME.message);
        fprintf('Continuing with sequential execution.\n');
    end
end

% Create all combinations of parameters
combinations = combvec(1:length(r_values), 1:length(variance_values))';
total_combinations = size(combinations, 1);

% Create parameter arrays that can be indexed inside parfor
r_idx_array = combinations(:, 1);
var_idx_array = combinations(:, 2);

% Progress reporting variables
progress_interval = max(1, round(total_combinations/20)); % Report progress ~20 times
fprintf('Total parameter combinations to process: %d\n', total_combinations);

% Create a temporary array to hold results (will be properly combined after parfor)
result_params = cell(total_combinations, 1);

% Use parfor for parallel execution
parfor i = 1:total_combinations
    % Extract parameter values for this iteration - use direct indexing 
    % which is compatible with parfor
    r_val = r_values(r_idx_array(i));
    var_val = variance_values(var_idx_array(i));
    
    % Create a local variable to store result
    current_result = [];
    
    try

        % Set random seed based on combination to ensure reproducibility
        % but different for each combination
        rng(); 
        
        %Generate lorenz data with varying initial conditions
        rng_condition = initial_conditions + randi([-5, 5], 1, 3); % Randomize initial conditions slightly
        [t,sol] = get_lorenz_data(rng_condition, tspan); % Call the function to generate data
        x_original = sol(:,1); % Use the first variable of the Lorenz system

        % Generate noisy data for this variance
        noise_std = sqrt(var_val); % Convert variance to standard deviation
        
        % Call the multi-run optimization function instead of the single-run version
        bestParams = findOptimalDenoiseParams_multi(x_original, noise_std, num_runs, r_val, dt, tspan, paramRanges);
        
        % Store the result in our temporary local variable
        current_result = bestParams;
        
        % Progress reporting (will be interleaved due to parallel execution)
        fprintf('Completed combination %d/%d: r = %d, variance = %g\n', ...
            i, total_combinations, r_val, var_val);
    catch ME
        fprintf('ERROR processing combination r=%d, variance=%g: %s\n', ...
            r_val, var_val, ME.message);
    end
    
    % Save the local result to our results array - this is compatible with parfor
    result_params{i} = current_result;
end

% Convert results to structured format
optimization_results = struct();

for i = 1:total_combinations
    % Get parameter values for this combination
    r_val = r_values(r_idx_array(i));
    var_val = variance_values(var_idx_array(i));
    bestParams = result_params{i};
    
    % Only process if we have results
    if ~isempty(bestParams)
        % Create valid field names
        r_field = sprintf('r%d', r_val);
        var_field = sprintf('var%g', var_val);
        var_field = strrep(var_field, '.', 'p')
        var_field = strrep(var_field, '-', 'neg');
        
        % Create structure fields as needed
        if ~isfield(optimization_results, r_field)
            optimization_results.(r_field) = struct();
        end
        
        % Store the results
        optimization_results.(r_field).(var_field) = bestParams;
    end
end

% Save the results structure to a .mat file
results_filename = 'optimal_denoise_params.mat';
try
    save(results_filename, 'optimization_results', 'paramRanges', 'r_values', 'variance_values', 'num_runs');
    fprintf('\n\nPrecomputation complete. Optimal parameters saved to %s\n', results_filename);
catch ME
    error('Failed to save results to %s. Error: %s', results_filename, ME.message);
end

% Display summary
fprintf('\nSummary of optimization results:\n');
fprintf('-------------------------------\n');
r_fields = fieldnames(optimization_results);
for i = 1:length(r_fields)
    r_field = r_fields{i};
    fprintf('For %s:\n', r_field);
    
    var_fields = fieldnames(optimization_results.(r_field));
    for j = 1:length(var_fields)
        var_field = var_fields{j};
        res = optimization_results.(r_field).(var_field);
        
        fprintf('  %s: Best method = %s (Error: %.6e)\n', ...
            var_field, res.best_method, res.best_error);
    end
end

fprintf('Total time taken: %.2f seconds\n', toc(clock_start));