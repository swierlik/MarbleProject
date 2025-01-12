% combined_generate_optimize_reconstruct.m

% Add this at the beginning for debugging
dbstop if error

% Parameters for Stochastic Lorenz System Generation
sigma = 10;
rho = 28;
beta = 8/3;

% Time span
dt = 0.001;
tfinal = 5; % Reduced from 50 to 5
tspan = dt:dt:tfinal; % Now length 5000

% Number of random systems to generate
n_systems = 1; % Set to 1 for debugging

% Adjust hankelRows
hankelRows = 50; % Reduced from 100 to 50

% Parameters for noise reduction optimization
ma_window_sizes = 3:2:21;
sg_orders = 2:4;
sg_framelens = 5:2:25;
wd_levels = 1:9;

% Preallocate results
best_params = zeros(n_systems, 3); % Columns: MA window size, SG order, WD level
system_errors = zeros(n_systems, 4); % Columns: Noisy, MA, SG, WD

% Loop for generating systems and finding best parameters
for sys_idx = 1:n_systems
    %% Generate Randomized Stochastic System
    % Randomized noise coefficients and initial conditions
    eta_x = 0.3 + 0.7 * rand; % Noise coefficient for x
    eta_y = 0.3 + 0.7 * rand;
    eta_z = 0.3 + 0.7 * rand;
    initial_conditions = -10 + 20 * rand(1, 3); % Random initial conditions

    % Additive stochastic noise
    dW_x = eta_x * sqrt(dt) * randn(size(tspan));
    dW_y = eta_y * sqrt(dt) * randn(size(tspan));
    dW_z = eta_z * sqrt(dt) * randn(size(tspan));

    % Initialize solution matrix (original deterministic data)
    deterministic_data = zeros(length(tspan), 3);
    deterministic_data(1, :) = initial_conditions;

    % Initialize solution matrix with noise
    noisy_data = zeros(length(tspan), 3);
    noisy_data(1, :) = initial_conditions;

    % Time integration for deterministic and noisy systems
    for i = 2:length(tspan)
        xyz = deterministic_data(i-1, :);
        deterministic = [sigma * (xyz(2) - xyz(1)); ...
                         xyz(1) * (rho - xyz(3)) - xyz(2); ...
                         xyz(1) * xyz(2) - beta * xyz(3)]';
        deterministic_data(i, :) = xyz + dt * deterministic;

        % Add stochastic noise for noisy data
        noisy_data(i, :) = deterministic_data(i, :) + [dW_x(i), dW_y(i), dW_z(i)];
    end

    %% Reconstruct Original Deterministic System
    try
        [~, A_det, B_det, xReg_det, t_reg_det] = getSystem(deterministic_data(:,1), hankelRows, 15, dt, tspan);

        % Check matrices for NaN or Inf
        if any(isnan(A_det(:))) || any(isinf(A_det(:))) || isempty(A_det)
            error('A_det contains NaN, Inf, or is empty');
        end
        if any(isnan(B_det(:))) || any(isinf(B_det(:))) || isempty(B_det)
            error('B_det contains NaN, Inf, or is empty');
        end

        sys_det = ss(A_det, B_det, eye(size(A_det, 1)), zeros(size(A_det,1), size(B_det,2)));
        [y_sim_det, ~] = lsim(sys_det, xReg_det(:, end), t_reg_det, xReg_det(1, 1:end-1));
    catch ME
        % Display detailed error information
        warning('Reconstruction failed for deterministic system at index %d: %s', sys_idx, ME.message);
        disp('Error identifier:');
        disp(ME.identifier);
        disp('Error stack:');
        for e = 1:length(ME.stack)
            disp(ME.stack(e));
        end
        continue;
    end

    %% Error Calculation and Optimization
    % Noisy reconstruction
    try
        [~, A_noisy, B_noisy, xReg_noisy, t_reg_noisy] = getSystem(noisy_data(:,1), hankelRows, 15, dt, tspan);

        % Check matrices for NaN or Inf
        if any(isnan(A_noisy(:))) || any(isinf(A_noisy(:))) || isempty(A_noisy)
            error('A_noisy contains NaN, Inf, or is empty');
        end
        if any(isnan(B_noisy(:))) || any(isinf(B_noisy(:))) || isempty(B_noisy)
            error('B_noisy contains NaN, Inf, or is empty');
        end

        sys_noisy = ss(A_noisy, B_noisy, eye(size(A_noisy, 1)), zeros(size(A_noisy,1), size(B_noisy,2)));
        [y_sim_noisy, ~] = lsim(sys_noisy, xReg_noisy(:, end), t_reg_noisy, xReg_noisy(1, 1:end-1));

        % Compute baseline error
        minLength = min(length(y_sim_det), length(y_sim_noisy));
        baseline_error = mean(sqrt(sum((y_sim_det(1:minLength,:) - y_sim_noisy(1:minLength,:)).^2, 2)));
    catch ME
        % If reconstruction fails, set baseline error to Inf
        baseline_error = Inf;
        warning('Reconstruction failed for noisy system at index %d: %s', sys_idx, ME.message);
        continue;
    end

    % Moving Average (MA)
    best_ma_error = Inf;
    best_ma_params = NaN;
    for j = 1:length(ma_window_sizes)
        window_size = ma_window_sizes(j);
        smoothed_data = movmean(noisy_data, window_size, 1); % Column-wise
        try
            [~, A_ma, B_ma, xReg_ma, t_reg_ma] = getSystem(smoothed_data(:,1), hankelRows, 15, dt, tspan);

            % Check matrices
            if any(isnan(A_ma(:))) || any(isinf(A_ma(:))) || isempty(A_ma)
                error('A_ma contains NaN, Inf, or is empty');
            end
            if any(isnan(B_ma(:))) || any(isinf(B_ma(:))) || isempty(B_ma)
                error('B_ma contains NaN, Inf, or is empty');
            end

            sys_ma = ss(A_ma, B_ma, eye(size(A_ma, 1)), zeros(size(A_ma,1), size(B_ma,2)));
            [y_sim_ma, ~] = lsim(sys_ma, xReg_ma(:, end), t_reg_ma, xReg_ma(1, 1:end-1));
            minLength = min(length(y_sim_det), length(y_sim_ma));
            error_ma = mean(sqrt(sum((y_sim_det(1:minLength,:) - y_sim_ma(1:minLength,:)).^2, 2)));
            if error_ma < best_ma_error
                best_ma_error = error_ma;
                best_ma_params = window_size;
            end
        catch
            continue;
        end
    end

    % Savitzky-Golay (SG)
    best_sg_error = Inf;
    best_sg_params = [NaN, NaN];
    for k = 1:length(sg_orders)
        for l = 1:length(sg_framelens)
            framelen = sg_framelens(l);
            if mod(framelen, 2) == 1 && framelen > sg_orders(k) + 2 && framelen <= size(noisy_data, 1)
                smoothed_data = sgolayfilt(noisy_data, sg_orders(k), framelen, [], 1); % Column-wise
                try
                    [~, A_sg, B_sg, xReg_sg, t_reg_sg] = getSystem(smoothed_data(:,1), hankelRows, 15, dt, tspan);

                    % Check matrices
                    if any(isnan(A_sg(:))) || any(isinf(A_sg(:))) || isempty(A_sg)
                        error('A_sg contains NaN, Inf, or is empty');
                    end
                    if any(isnan(B_sg(:))) || any(isinf(B_sg(:))) || isempty(B_sg)
                        error('B_sg contains NaN, Inf, or is empty');
                    end

                    sys_sg = ss(A_sg, B_sg, eye(size(A_sg, 1)), zeros(size(A_sg,1), size(B_sg,2)));
                    [y_sim_sg, ~] = lsim(sys_sg, xReg_sg(:, end), t_reg_sg, xReg_sg(1, 1:end-1));
                    minLength = min(length(y_sim_det), length(y_sim_sg));
                    error_sg = mean(sqrt(sum((y_sim_det(1:minLength,:) - y_sim_sg(1:minLength,:)).^2, 2)));
                    if error_sg < best_sg_error
                        best_sg_error = error_sg;
                        best_sg_params = [sg_orders(k), framelen];
                    end
                catch
                    continue;
                end
            end
        end
    end

    % Wavelet Denoising (WD)
    best_wd_error = Inf;
    best_wd_params = NaN;
    for m = 1:length(wd_levels)
        level = wd_levels(m);
        smoothed_data = zeros(size(noisy_data)); % Preallocate
        for col = 1:size(noisy_data, 2)
            smoothed_data(:, col) = wdenoise(noisy_data(:, col), level);
        end
        try
            [~, A_wd, B_wd, xReg_wd, t_reg_wd] = getSystem(smoothed_data(:,1), hankelRows, 15, dt, tspan);

            % Check matrices
            if any(isnan(A_wd(:))) || any(isinf(A_wd(:))) || isempty(A_wd)
                error('A_wd contains NaN, Inf, or is empty');
            end
            if any(isnan(B_wd(:))) || any(isinf(B_wd(:))) || isempty(B_wd)
                error('B_wd contains NaN, Inf, or is empty');
            end

            sys_wd = ss(A_wd, B_wd, eye(size(A_wd, 1)), zeros(size(A_wd,1), size(B_wd,2)));
            [y_sim_wd, ~] = lsim(sys_wd, xReg_wd(:, end), t_reg_wd, xReg_wd(1, 1:end-1));
            minLength = min(length(y_sim_det), length(y_sim_wd));
            error_wd = mean(sqrt(sum((y_sim_det(1:minLength,:) - y_sim_wd(1:minLength,:)).^2, 2)));
            if error_wd < best_wd_error
                best_wd_error = error_wd;
                best_wd_params = level;
            end
        catch
            continue;
        end
    end

    % Save best parameters and errors
    best_params(sys_idx, :) = [best_ma_params, best_sg_params(1), best_wd_params];
    system_errors(sys_idx, :) = [baseline_error, best_ma_error, best_sg_error, best_wd_error];
end

%% Aggregate Results
% Remove any NaN or Inf values from results
valid_indices = all(~isinf(system_errors) & ~isnan(system_errors), 2);
best_params = best_params(valid_indices, :);
system_errors = system_errors(valid_indices, :);

% Find the most common best parameters across systems
if ~isempty(best_params)
    final_best_ma = mode(best_params(:, 1));
    final_best_sg_order = mode(best_params(:, 2));
    final_best_wd_level = mode(best_params(:, 3));
else
    final_best_ma = NaN;
    final_best_sg_order = NaN;
    final_best_wd_level = NaN;
end

% Print universal best parameters
fprintf('Universal Best Moving Average Window Size: %d\n', final_best_ma);
fprintf('Universal Best Savitzky-Golay Order: %d\n', final_best_sg_order);
fprintf('Universal Best Wavelet Denoising Level: %d\n', final_best_wd_level);

% Plot system errors for visualization
if ~isempty(system_errors)
    figure;
    bar(mean(system_errors, 1));
    xticks(1:4);
    xticklabels({'Noisy', 'MA', 'SG', 'WD'});
    title('Average Reconstruction Errors Across Systems');
    ylabel('Mean Error');
    grid on;
else
    disp('No valid system errors to display.');
end
