% noise_optimization_multiple_realizations.m

% Load the data
load('lorenzData.mat') % Contains 'sol', 't', 'dt'
load('systemData.mat') % Contains 'V_x', 'V_y', 'A_x', 'B_x', 'A_y', 'B_y', 'xReg', 'yReg', 'r', 'tspan', 'dt')

% Extract x, y, z from sol
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Define L (exclude initial and final transients)
L = 300:length(xReg)-300;
forcing_vector = xReg(L, r);

% Determine the number of input channels and define the time vector
n_inputs = size(B_x, 2);
time_vector = dt * (L - 1);

% Set the seed for reproducibility
rng(117);

% Number of noisy realizations
n_realizations = 100;

% Preallocate error tracking arrays
errors_noisy = zeros(n_realizations, length(time_vector));
errors_ma = zeros(n_realizations, length(time_vector));
errors_sg = zeros(n_realizations, length(time_vector));
errors_wd = zeros(n_realizations, length(time_vector));

% Initialize best parameters for each method
best_ma_window_size = 5;
best_sg_order = 3;
best_sg_framelen = 7;
best_wd_level = 3;
min_mean_error_ma = Inf;
min_mean_error_sg = Inf;
min_mean_error_wd = Inf;

% Define ranges for parameter search
ma_window_sizes = 3:2:21;
sg_orders = 2:4;
sg_framelens = 5:2:25;
wd_levels = 1:9;

% Loop over all noisy realizations
for i = 1:n_realizations
    % Generate a new noisy vector
    noise_level = 0.1 * std(forcing_vector);
    noisy_forcing_vector = forcing_vector + noise_level * randn(size(forcing_vector));
    noisy_forcing_matrix = repmat(noisy_forcing_vector, 1, n_inputs);
    [y_sim_x_noisy, ~] = lsim(sys_x, noisy_forcing_matrix, time_vector, xReg(L(1), 1:r-1));
    error_noisy = sqrt(sum((y_sim_x_noisy - y_sim_x).^2, 2));
    errors_noisy(i, :) = error_noisy;

    % ------------------ Optimization for Moving Average Filter ------------------
    for window_size = ma_window_sizes
        smoothed_forcing_vector_ma = movmean(noisy_forcing_vector, window_size);
        smoothed_forcing_matrix_ma = repmat(smoothed_forcing_vector_ma, 1, n_inputs);
        [y_sim_x_ma, ~] = lsim(sys_x, smoothed_forcing_matrix_ma, time_vector, xReg(L(1), 1:r-1));
        error_ma = sqrt(sum((y_sim_x_ma - y_sim_x).^2, 2));

        % Update best parameters if average error improves
        mean_error_ma = mean(error_ma);
        if mean_error_ma < min_mean_error_ma
            min_mean_error_ma = mean_error_ma;
            best_ma_window_size = window_size;
        end
    end

    % Store the error for the best Moving Average parameters
    smoothed_forcing_vector_ma_best = movmean(noisy_forcing_vector, best_ma_window_size);
    smoothed_forcing_matrix_ma_best = repmat(smoothed_forcing_vector_ma_best, 1, n_inputs);
    [y_sim_x_ma_best, ~] = lsim(sys_x, smoothed_forcing_matrix_ma_best, time_vector, xReg(L(1), 1:r-1));
    errors_ma(i, :) = sqrt(sum((y_sim_x_ma_best - y_sim_x).^2, 2));

    % ------------------ Optimization for Savitzky-Golay Filter ------------------
    for order = sg_orders
        for framelen = sg_framelens
            if mod(framelen, 2) == 1
                smoothed_forcing_vector_sg = sgolayfilt(noisy_forcing_vector, order, framelen);
                smoothed_forcing_matrix_sg = repmat(smoothed_forcing_vector_sg, 1, n_inputs);
                [y_sim_x_sg, ~] = lsim(sys_x, smoothed_forcing_matrix_sg, time_vector, xReg(L(1), 1:r-1));
                error_sg = sqrt(sum((y_sim_x_sg - y_sim_x).^2, 2));

                mean_error_sg = mean(error_sg);
                if mean_error_sg < min_mean_error_sg
                    min_mean_error_sg = mean_error_sg;
                    best_sg_order = order;
                    best_sg_framelen = framelen;
                end
            end
        end
    end

    % Store the error for the best Savitzky-Golay parameters
    smoothed_forcing_vector_sg_best = sgolayfilt(noisy_forcing_vector, best_sg_order, best_sg_framelen);
    smoothed_forcing_matrix_sg_best = repmat(smoothed_forcing_vector_sg_best, 1, n_inputs);
    [y_sim_x_sg_best, ~] = lsim(sys_x, smoothed_forcing_matrix_sg_best, time_vector, xReg(L(1), 1:r-1));
    errors_sg(i, :) = sqrt(sum((y_sim_x_sg_best - y_sim_x).^2, 2));

    % ------------------ Optimization for Wavelet Denoising ------------------
    for level = wd_levels
        smoothed_forcing_vector_wd = wdenoise(noisy_forcing_vector, level);
        smoothed_forcing_matrix_wd = repmat(smoothed_forcing_vector_wd, 1, n_inputs);
        [y_sim_x_wd, ~] = lsim(sys_x, smoothed_forcing_matrix_wd, time_vector, xReg(L(1), 1:r-1));
        error_wd = sqrt(sum((y_sim_x_wd - y_sim_x).^2, 2));

        mean_error_wd = mean(error_wd);
        if mean_error_wd < min_mean_error_wd
            min_mean_error_wd = mean_error_wd;
            best_wd_level = level;
        end
    end

    % Store the error for the best Wavelet Denoising parameters
    smoothed_forcing_vector_wd_best = wdenoise(noisy_forcing_vector, best_wd_level);
    smoothed_forcing_matrix_wd_best = repmat(smoothed_forcing_vector_wd_best, 1, n_inputs);
    [y_sim_x_wd_best, ~] = lsim(sys_x, smoothed_forcing_matrix_wd_best, time_vector, xReg(L(1), 1:r-1));
    errors_wd(i, :) = sqrt(sum((y_sim_x_wd_best - y_sim_x).^2, 2));
end

% ------------------ Final Time-Series Error Plot ------------------

% Compute the mean error across all realizations
mean_error_noisy = mean(errors_noisy, 1);
mean_error_ma = mean(errors_ma, 1);
mean_error_sg = mean(errors_sg, 1);
mean_error_wd = mean(errors_wd, 1);

figure;
plot(time_vector, mean_error_noisy, 'r', 'LineWidth', 1);
hold on;
plot(time_vector, mean_error_ma, 'g', 'LineWidth', 1);
plot(time_vector, mean_error_sg, 'm', 'LineWidth', 1);
plot(time_vector, mean_error_wd, 'k', 'LineWidth', 1);
legend('Noisy', 'MA Optimized', 'SG Optimized', 'WD Optimized');
title('Mean Error Over Time (Averaged Across 100 Realizations)');
xlabel('Time');
ylabel('Mean Error');
grid on;

% ------------------ End of Script ------------------
