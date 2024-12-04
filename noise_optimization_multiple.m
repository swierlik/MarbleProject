% noise_optimization_multiple_stochastic.m

% Load the data
load('Data/lorenzData.mat'); % Base deterministic Lorenz data
load('Data/lorenzDataStochastic.mat'); % Stochastic Lorenz data
load('Data/systemData.mat'); % Contains 'V_x', 'V_y', 'A_x', 'B_x', 'A_y', 'B_y', 'xReg', 'yReg', 'r', 'tspan', 'dt'

% Extract x, y, z from base deterministic solution
x_base = sol(:, 1);
y_base = sol(:, 2);
z_base = sol(:, 3);

% Extract x, y, z from stochastic solution
x_stochastic = sol2(:, 1);
y_stochastic = sol2(:, 2);
z_stochastic = sol2(:, 3);

% Ensure the lengths of stochastic and deterministic data match
min_length = min(length(x_base), length(x_stochastic));
x_base = x_base(1:min_length);
y_base = y_base(1:min_length);
z_base = z_base(1:min_length);

x_stochastic = x_stochastic(1:min_length);
y_stochastic = y_stochastic(1:min_length);
z_stochastic = z_stochastic(1:min_length);

% Combine into matrices
deterministic_data = [x_base, y_base, z_base];
stochastic_data = [x_stochastic, y_stochastic, z_stochastic];

% Define L (exclude initial and final transients)
L = 300:(min_length - 300);

% Time vector
time_vector = dt * (L - 1);

% Number of noisy realizations
n_realizations = size(stochastic_data, 1);

% Define ranges for parameter search
ma_window_sizes = 3:2:21;
sg_orders = 2:4;
sg_framelens = 5:2:25;
wd_levels = 1:9;

% Preallocate for average errors
avg_error_noisy = 0;
avg_error_ma = zeros(1, length(ma_window_sizes));
avg_error_sg = zeros(length(sg_orders), length(sg_framelens));
avg_error_wd = zeros(1, length(wd_levels));

% Debugging output: Check data sizes
fprintf('Number of realizations: %d\n', n_realizations);
fprintf('Data size (stochastic): [%d, %d]\n', size(stochastic_data));
fprintf('Data size (deterministic): [%d, %d]\n', size(deterministic_data));

% Parallelized loop for realizations
parfor i = 1:size(stochastic_data, 1)
    % Temporary variables for each worker
    temp_error_noisy = 0;
    temp_error_ma = zeros(1, length(ma_window_sizes));
    temp_error_sg = zeros(length(sg_orders), length(sg_framelens));
    temp_error_wd = zeros(1, length(wd_levels));

    % Extract current noisy realization
    noisy_data = stochastic_data(i, :);

    % Compute noisy error (baseline)
    temp_error_noisy = sqrt(mean((noisy_data - deterministic_data(i, :)).^2));

    % ------------------ Optimization for Moving Average Filter ------------------
    for j = 1:length(ma_window_sizes)
        window_size = ma_window_sizes(j);
        smoothed_data_ma = zeros(size(noisy_data)); % Preallocate for smoothed data
        for col = 1:size(noisy_data, 2) % Apply column-wise
            smoothed_data_ma(:, col) = movmean(noisy_data(:, col), window_size);
        end
        error_ma = sqrt(mean((smoothed_data_ma - deterministic_data(i, :)).^2));
        temp_error_ma(j) = temp_error_ma(j) + error_ma;
    end

    % ------------------ Optimization for Savitzky-Golay Filter ------------------
    for k = 1:length(sg_orders)
        for l = 1:length(sg_framelens)
            framelen = sg_framelens(l);
            if mod(framelen, 2) == 1 && framelen <= length(noisy_data) % Ensure valid frame length
                order = sg_orders(k);
                smoothed_data_sg = zeros(size(noisy_data)); % Preallocate for smoothed data
                for col = 1:size(noisy_data, 2) % Apply column-wise
                    smoothed_data_sg(:, col) = sgolayfilt(noisy_data(:, col), order, framelen);
                end
                error_sg = sqrt(mean((smoothed_data_sg - deterministic_data(i, :)).^2));
                temp_error_sg(k, l) = temp_error_sg(k, l) + error_sg;
            end
        end
    end

    % ------------------ Optimization for Wavelet Denoising ------------------
    for m = 1:length(wd_levels)
        level = wd_levels(m);
        smoothed_data_wd = zeros(size(noisy_data)); % Preallocate for smoothed data
        for col = 1:size(noisy_data, 2) % Apply column-wise
            if length(noisy_data(:, col)) > 1 % Ensure valid input length
                smoothed_data_wd(:, col) = wdenoise(noisy_data(:, col), level);
            else
                smoothed_data_wd(:, col) = noisy_data(:, col); % Leave unchanged if too short
            end
        end
        error_wd = sqrt(mean((smoothed_data_wd - deterministic_data(i, :)).^2));
        temp_error_wd(m) = temp_error_wd(m) + error_wd;
    end

    % Aggregate temporary results
    avg_error_noisy = avg_error_noisy + temp_error_noisy;
    avg_error_ma = avg_error_ma + temp_error_ma;
    avg_error_sg = avg_error_sg + temp_error_sg;
    avg_error_wd = avg_error_wd + temp_error_wd;
end

% ------------------ Normalize Average Errors ------------------
avg_error_noisy = avg_error_noisy / n_realizations;
avg_error_ma = avg_error_ma / n_realizations;
avg_error_sg = avg_error_sg / n_realizations;
avg_error_wd = avg_error_wd / n_realizations;

% ------------------ Determine Best Parameters ------------------

% Find the best parameters for each method
[~, best_ma_idx] = min(avg_error_ma);
best_ma_window_size = ma_window_sizes(best_ma_idx);

[min_sg_error, sg_idx] = min(avg_error_sg(:));
[best_sg_order_idx, best_sg_framelen_idx] = ind2sub(size(avg_error_sg), sg_idx);
best_sg_order = sg_orders(best_sg_order_idx);
best_sg_framelen = sg_framelens(best_sg_framelen_idx);

[~, best_wd_idx] = min(avg_error_wd);
best_wd_level = wd_levels(best_wd_idx);

% ------------------ Print Best Parameters ------------------
fprintf('Best Moving Average Window Size: %d\n', best_ma_window_size);
fprintf('Best Savitzky-Golay Order: %d, Frame Length: %d\n', best_sg_order, best_sg_framelen);
fprintf('Best Wavelet Denoising Level: %d\n', best_wd_level);
