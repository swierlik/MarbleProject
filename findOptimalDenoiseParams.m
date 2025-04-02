% findOptimalDenoiseParams.m
function bestParams = findOptimalDenoiseParams(x_original, x_noisy, r, dt, tspan, paramRanges)
    % Finds optimal denoising parameters for Moving Average, Savitzky-Golay,
    % and Wavelet methods for a given noise level and HAVOK rank 'r'.
    %
    % Inputs:
    %   x_original  - Original clean signal
    %   x_noisy     - Noisy signal
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
    
        fprintf('--- Optimizing for r=%d ---\n', r);
    
        % Initialize results containers for this run
        best_ma_error = Inf;
        best_ma_window = paramRanges.window_sizes(1); % Default fallback
    
        best_sg_error = Inf;
        best_sg_order = paramRanges.sg_orders(1);     % Default fallback
        best_sg_framelen = paramRanges.sg_framelens(1); % Default fallback Adjusted below
    
        best_wd_error = Inf;
        best_wd_level = paramRanges.wd_levels(1);     % Default fallback
        best_wd_wavelet = paramRanges.wd_wavelets{1}; % Default fallback
    
        % Adjust SG framelen default if needed
        if best_sg_framelen <= best_sg_order
            valid_framerates = paramRanges.sg_framelens(paramRanges.sg_framelens > best_sg_order & mod(paramRanges.sg_framelens, 2) ~= 0);
            if ~isempty(valid_framerates)
                best_sg_framelen = valid_framerates(1);
            else
                 warning('No valid SG framelen found for default order %d.', best_sg_order);
                 % Handle error appropriately, maybe skip SG or set error to Inf
                 best_sg_error = Inf; % Indicate SG couldn't run
            end
        end
    
    
        % -------- Optimize Moving Average Parameters --------
        fprintf('Optimizing Moving Average...\n');
        for window_size = paramRanges.window_sizes
            fprintf('  Testing MA window size: %d\n', window_size);
            try
                x_movmean = movmean(x_noisy, window_size);
                error_val = compute_reconstruction_error(x_movmean, x_original, r, dt, tspan);
                if error_val < best_ma_error
                    best_ma_error = error_val;
                    best_ma_window = window_size;
                    fprintf('  New best MA window: %d (Error: %.6e)\n', window_size, error_val);
                end
            catch ME
                 fprintf('  Error during MA calculation for window %d: %s\n', window_size, ME.message);
            end
        end
    
        % -------- Optimize Savitzky-Golay Parameters --------
        fprintf('Optimizing Savitzky-Golay...\n');
        for order = paramRanges.sg_orders
            valid_framerates = paramRanges.sg_framelens(paramRanges.sg_framelens > order & mod(paramRanges.sg_framelens, 2) ~= 0);
            if isempty(valid_framerates)
                 fprintf('  No valid SG framelen found for order %d. Skipping order.\n', order);
                 continue;
            end
            for framelen = valid_framerates
                fprintf('  Testing SG order: %d, framelen: %d\n', order, framelen);
                 try
                    x_sg = sgolayfilt(x_noisy, order, framelen);
                    error_val = compute_reconstruction_error(x_sg, x_original, r, dt, tspan);
                    if error_val < best_sg_error
                        best_sg_error = error_val;
                        best_sg_order = order;
                        best_sg_framelen = framelen;
                        fprintf('  New best SG: order=%d, framelen=%d (Error: %.6e)\n', order, framelen, error_val);
                    end
                 catch ME
                     fprintf('  Error during SG calculation for order %d, framelen %d: %s\n', order, framelen, ME.message);
                 end
            end
        end
    
        % -------- Optimize Wavelet Denoising Parameters --------
        fprintf('Optimizing Wavelet Denoising...\n');
        for level = paramRanges.wd_levels
            for w = 1:length(paramRanges.wd_wavelets)
                wavelet = paramRanges.wd_wavelets{w};
                fprintf('  Testing WD level: %d, wavelet: %s\n', level, wavelet);
                try
                    x_wd = wdenoise(x_noisy, level, 'Wavelet', wavelet);
                    error_val = compute_reconstruction_error(x_wd, x_original, r, dt, tspan);
                    if error_val < best_wd_error
                        best_wd_error = error_val;
                        best_wd_level = level;
                        best_wd_wavelet = wavelet;
                        fprintf('  New best WD: level=%d, wavelet=%s (Error: %.6e)\n', level, wavelet, error_val);
                    end
                catch ME
                    fprintf('  Error during WD calculation for level %d, wavelet %s: %s\n', level, wavelet, ME.message);
                    % Consider skipping this wavelet if it consistently fails
                end
            end
        end
    
        % Store best parameters found
        bestParams.ma_window_size = best_ma_window;
        bestParams.sg_order = best_sg_order;
        bestParams.sg_framelen = best_sg_framelen;
        bestParams.wd_level = best_wd_level;
        bestParams.wd_wavelet = best_wd_wavelet;
    
        fprintf('--- Optimization complete for r=%d ---\n', r);
    
    end % End of findOptimalDenoiseParams function
    
    
    % --- NESTED HELPER FUNCTION ---
    function error_val = compute_reconstruction_error(x_denoised, x_original_nested, r_nested, dt_nested, tspan_nested)
    % Computes the mean squared error of the HAVOK reconstruction.
    % Uses variables from the parent function's scope.
    
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