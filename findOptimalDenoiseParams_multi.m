% findOptimalDenoiseParams_multi.m
function bestParams = findOptimalDenoiseParams_multi(x_original, noise_std, num_runs, r, dt, tspan, paramRanges)
    % Finds optimal denoising parameters by averaging results over multiple noise realizations
    %
    % Inputs:
    %   x_original  - Original clean signal
    %   noise_std   - Standard deviation of noise to add
    %   num_runs    - Number of noise seeds/realizations to average over
    %   r           - HAVOK rank parameter
    %   dt          - Time step
    %   tspan       - Time vector for simulation
    %   paramRanges - Struct containing parameter ranges to test:
    %                 .window_sizes
    %                 .sg_orders
    %                 .sg_framelens
    %                 .wd_levels
    %                 .wd_wavelets
    %
    % Output:
    %   bestParams  - Struct containing the optimal parameters found:
    %                 .ma_window_size
    %                 .sg_order
    %                 .sg_framelen
    %                 .wd_level
    %                 .wd_wavelet
    %                 .average_errors (struct with errors for each method)
    
    fprintf('=== Finding optimal parameters over %d noise realizations (r=%d) ===\n', num_runs, r);

    rng(); % Reset random seed for each function call
    
    % Initialize error accumulators for all parameter combinations
    ma_errors = zeros(length(paramRanges.window_sizes), num_runs);
    sg_errors = zeros(length(paramRanges.sg_orders), length(paramRanges.sg_framelens), num_runs);
    wd_errors = zeros(length(paramRanges.wd_levels), length(paramRanges.wd_wavelets), num_runs);
    
    % Initialize parameter indices for tracking
    ma_params_idx = zeros(length(paramRanges.window_sizes), 1);
    for i = 1:length(paramRanges.window_sizes)
        ma_params_idx(i) = paramRanges.window_sizes(i);
    end
    
    sg_order_idx = zeros(length(paramRanges.sg_orders), 1);
    for i = 1:length(paramRanges.sg_orders)
        sg_order_idx(i) = paramRanges.sg_orders(i);
    end
    
    sg_frame_idx = zeros(length(paramRanges.sg_framelens), 1);
    for i = 1:length(paramRanges.sg_framelens)
        sg_frame_idx(i) = paramRanges.sg_framelens(i);
    end
    
    % Run multiple trials with different noise seeds
    for run = 1:num_runs
        fprintf('\n--- Run %d/%d ---\n', run, num_runs);
        
        % Generate noisy signal for this run
        noise = noise_std * randn(size(x_original));
        x_noisy = x_original + noise;
        
        % -------- Test Moving Average Parameters --------
        fprintf('Testing Moving Average parameters...\n');
        for i = 1:length(paramRanges.window_sizes)
            window_size = paramRanges.window_sizes(i);
            fprintf('  Testing MA window size: %d\n', window_size);
            try
                x_movmean = movmean(x_noisy, window_size);
                error_val = compute_reconstruction_error(x_movmean, x_original, r, dt, tspan);
                ma_errors(i, run) = error_val;
                fprintf('  MA window %d error: %.6e\n', window_size, error_val);
            catch ME
                fprintf('  Error during MA calculation for window %d: %s\n', window_size, ME.message);
                ma_errors(i, run) = Inf;
            end
        end
        
        % -------- Test Savitzky-Golay Parameters --------
        fprintf('Testing Savitzky-Golay parameters...\n');
        for i = 1:length(paramRanges.sg_orders)
            order = paramRanges.sg_orders(i);
            for j = 1:length(paramRanges.sg_framelens)
                framelen = paramRanges.sg_framelens(j);
                
                % SG filter requires that framelen > order and framelen is odd
                if framelen <= order || mod(framelen, 2) == 0
                    sg_errors(i, j, run) = Inf;
                    continue;
                end
                
                fprintf('  Testing SG order: %d, framelen: %d\n', order, framelen);
                try
                    x_sg = sgolayfilt(x_noisy, order, framelen);
                    error_val = compute_reconstruction_error(x_sg, x_original, r, dt, tspan);
                    sg_errors(i, j, run) = error_val;
                    fprintf('  SG order %d, framelen %d error: %.6e\n', order, framelen, error_val);
                catch ME
                    fprintf('  Error during SG calculation for order %d, framelen %d: %s\n', order, framelen, ME.message);
                    sg_errors(i, j, run) = Inf;
                end
            end
        end
        
        % -------- Test Wavelet Denoising Parameters --------
        fprintf('Testing Wavelet Denoising parameters...\n');
        for i = 1:length(paramRanges.wd_levels)
            level = paramRanges.wd_levels(i);
            for j = 1:length(paramRanges.wd_wavelets)
                wavelet = paramRanges.wd_wavelets{j};
                fprintf('  Testing WD level: %d, wavelet: %s\n', level, wavelet);
                try
                    x_wd = wdenoise(x_noisy, level, 'Wavelet', wavelet);
                    error_val = compute_reconstruction_error(x_wd, x_original, r, dt, tspan);
                    wd_errors(i, j, run) = error_val;
                    fprintf('  WD level %d, wavelet %s error: %.6e\n', level, wavelet, error_val);
                catch ME
                    fprintf('  Error during WD calculation for level %d, wavelet %s: %s\n', level, wavelet, ME.message);
                    wd_errors(i, j, run) = Inf;
                end
            end
        end
    end % End of runs loop
    
    % Calculate average errors across all runs
    ma_avg_errors = mean(ma_errors, 2);
    sg_avg_errors = mean(sg_errors, 3);
    wd_avg_errors = mean(wd_errors, 3);
    
    % Find best parameters based on average errors
    [min_ma_error, min_ma_idx] = min(ma_avg_errors);
    best_ma_window = paramRanges.window_sizes(min_ma_idx);
    
    % For SG, need to find the minimum in a 2D matrix
    [min_sg_error, linear_idx] = min(sg_avg_errors(:));
    [i, j] = ind2sub(size(sg_avg_errors), linear_idx);
    best_sg_order = paramRanges.sg_orders(i);
    best_sg_framelen = paramRanges.sg_framelens(j);
    
    % For WD, need to find the minimum in a 2D matrix
    [min_wd_error, linear_idx] = min(wd_avg_errors(:));
    [i, j] = ind2sub(size(wd_avg_errors), linear_idx);
    best_wd_level = paramRanges.wd_levels(i);
    best_wd_wavelet = paramRanges.wd_wavelets{j};
    
    % Store best parameters found and their average errors
    bestParams.ma_window_size = best_ma_window;
    bestParams.sg_order = best_sg_order;
    bestParams.sg_framelen = best_sg_framelen;
    bestParams.wd_level = best_wd_level;
    bestParams.wd_wavelet = best_wd_wavelet;
    
    % Store the error values
    bestParams.average_errors.ma = min_ma_error;
    bestParams.average_errors.sg = min_sg_error;
    bestParams.average_errors.wd = min_wd_error;
    
    % Create detailed report
    fprintf('\n=== OPTIMIZATION RESULTS (Averaged over %d runs) ===\n', num_runs);
    fprintf('Best Moving Average: window_size = %d (avg error: %.6e)\n', best_ma_window, min_ma_error);
    fprintf('Best Savitzky-Golay: order = %d, framelen = %d (avg error: %.6e)\n', best_sg_order, best_sg_framelen, min_sg_error);
    fprintf('Best Wavelet Denoising: level = %d, wavelet = %s (avg error: %.6e)\n', best_wd_level, best_wd_wavelet, min_wd_error);
    
    % Determine overall best method
    [best_error, best_method_idx] = min([min_ma_error, min_sg_error, min_wd_error]);
    methods = {'Moving Average', 'Savitzky-Golay', 'Wavelet Denoising'};
    bestParams.best_method = methods{best_method_idx};
    bestParams.best_error = best_error;
    
    fprintf('\nOVERALL BEST METHOD: %s (avg error: %.6e)\n', bestParams.best_method, best_error);
    fprintf('=== Optimization complete ===\n\n');
end % End of findOptimalDenoiseParams_multi function


% --- HELPER FUNCTION ---
function error_val = compute_reconstruction_error(x_denoised, x_original_nested, r_nested, dt_nested, tspan_nested)
    % Computes the mean squared error of the HAVOK reconstruction.
    % This function is identical to the original one
    
    error_val = Inf; % Default error if calculation fails
    
    try
        % Get HAVOK system for the denoised data
        % Assuming getSystem & dehankelize are available on the path
        [V, A, B, xReg, U, E] = getSystem(x_denoised, 100, r_nested, dt_nested, tspan_nested);
        
        % Check if getSystem produced valid outputs
        if isempty(A) || isempty(B) || isempty(xReg) || isempty(U) || isempty(E) || size(A,1) ~= r_nested-1
            warning('getSystem did not return valid matrices for r=%d. Skipping error calculation.', r_nested);
            return;
        end
        
        % Simulate the system
        L = 1:min(length(tspan_nested), size(xReg, 1));
        if isempty(L) || r_nested <= 1
            warning('Invalid simulation range or r value. Skipping error calculation.');
            return;
        end
        
        % Ensure U and E match Koopman mode rank r-1
        if size(U, 2) < r_nested-1 || size(E, 1) < r_nested-1 || size(E, 2) < r_nested-1
            warning('U or E matrix dimensions (%dx%d, %dx%d) are smaller than required r-1 (%d). Skipping error calculation.', ...
                size(U,1), size(U,2), size(E,1), size(E,2), r_nested-1);
            return;
        end
        U = U(:, 1:r_nested-1);
        E = E(1:r_nested-1, 1:r_nested-1);
        
        % System matrices check
        if size(A,1) ~= size(B,1) || size(A,1) ~= (r_nested-1) || size(B,2) ~= 1
            warning('System matrix dimensions mismatch (A:%dx%d, B:%dx%d, r-1=%d). Skipping simulation.', ...
                size(A,1), size(A,2), size(B,1), size(B,2), r_nested-1);
            return; % Cannot simulate
        end
        
        sys = ss(A, B, eye(r_nested-1), 0*B);
        
        % Ensure initial condition and input signal dimensions match
        initial_condition = xReg(1, 1:r_nested-1).'; % Ensure column vector
        input_signal = xReg(L, r_nested);
        time_vector = dt_nested*(L-1);
        
        if size(initial_condition, 1) ~= (r_nested-1)
            warning('Initial condition size mismatch. Skipping simulation.');
            return;
        end
        if isempty(input_signal)
            warning('Input signal is empty. Skipping simulation.');
            return;
        end
        
        
        [y_sim, ~] = lsim(sys, input_signal, time_vector, initial_condition);
        
        if isempty(y_sim) || size(y_sim, 2) ~= (r_nested-1)
            warning('Simulation failed or returned unexpected dimensions. Skipping error calculation.');
            return;
        end
        
        % Reconstruct
        x_reconstructed = U * E * y_sim.'; % y_sim is N x (r-1) -> y_sim.' is (r-1) x N
        x_reconstructed = dehankelize(x_reconstructed); % Assuming dehankelize returns a column vector
        
        % Ensure dimensions match for error calculation
        len_recon = length(x_reconstructed);
        len_orig = length(x_original_nested);
        len_L = length(L);
        
        % Determine the valid comparison length
        compare_len = min([len_recon, len_orig, len_L]);
        
        if compare_len < 1
            warning('Cannot compare signals due to length mismatch or zero length.');
            return;
        end
        
        % Compute error relative to original data over the valid range
        error_val = mean((x_original_nested(1:compare_len) - x_reconstructed(1:compare_len)).^2);
        
    catch ME
        warning('Error during reconstruction/simulation for r=%d: %s', r_nested, ME.message);
        fprintf('Error occurred in file %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
        error_val = Inf; % Assign high error if anything fails
    end
end