% noise_optimization.m

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

%Setting the seed for reproductability
rng(117);

% Simulate original system (without noise) for comparison
sys_x = ss(A_x, B_x, eye(r-1), zeros(r-1, n_inputs));
forcing_matrix = repmat(forcing_vector, 1, n_inputs);
[y_sim_x, t_sim_x] = lsim(sys_x, forcing_matrix, time_vector, xReg(L(1), 1:r-1));

% Add Gaussian noise
noise_level = 0.1 * std(forcing_vector);
noisy_forcing_vector = forcing_vector + noise_level * randn(size(forcing_vector));
noisy_forcing_matrix = repmat(noisy_forcing_vector, 1, n_inputs);
[y_sim_x_noisy, ~] = lsim(sys_x, noisy_forcing_matrix, time_vector, xReg(L(1), 1:r-1));
error_noisy = sqrt(sum((y_sim_x_noisy - y_sim_x).^2, 2));

% ------------------ Optimization for Moving Average Filter ------------------

initial_window_size = 5;
best_ma_window_size = initial_window_size;
min_error_ma = Inf;

for window_size = 3:2:21
    smoothed_forcing_vector_ma = movmean(noisy_forcing_vector, window_size);
    smoothed_forcing_matrix_ma = repmat(smoothed_forcing_vector_ma, 1, n_inputs);
    [y_sim_x_ma, ~] = lsim(sys_x, smoothed_forcing_matrix_ma, time_vector, xReg(L(1), 1:r-1));
    error_ma = sqrt(sum((y_sim_x_ma - y_sim_x).^2, 2));
    mean_error_ma = mean(error_ma);

    if mean_error_ma < min_error_ma
        min_error_ma = mean_error_ma;
        best_ma_window_size = window_size;
    end
end

% Initial and optimized vectors for Moving Average
smoothed_forcing_vector_ma_initial = movmean(noisy_forcing_vector, initial_window_size);
smoothed_forcing_vector_ma_best = movmean(noisy_forcing_vector, best_ma_window_size);

% ------------------ Optimization for Savitzky-Golay Filter ------------------

initial_order = 3;
initial_framelen = 7;
best_sg_order = initial_order;
best_sg_framelen = initial_framelen;
min_error_sg = Inf;

for order = 2:4
    for framelen = max(5, order+2):2:50
        if mod(framelen, 2) == 1
            smoothed_forcing_vector_sg = sgolayfilt(noisy_forcing_vector, order, framelen);
            smoothed_forcing_matrix_sg = repmat(smoothed_forcing_vector_sg, 1, n_inputs);
            [y_sim_x_sg, ~] = lsim(sys_x, smoothed_forcing_matrix_sg, time_vector, xReg(L(1), 1:r-1));
            error_sg = sqrt(sum((y_sim_x_sg - y_sim_x).^2, 2));
            mean_error_sg = mean(error_sg);

            if mean_error_sg < min_error_sg
                min_error_sg = mean_error_sg;
                best_sg_order = order;
                best_sg_framelen = framelen;
            end
        end
    end
end

smoothed_forcing_vector_sg_initial = sgolayfilt(noisy_forcing_vector, initial_order, initial_framelen);
smoothed_forcing_vector_sg_best = sgolayfilt(noisy_forcing_vector, best_sg_order, best_sg_framelen);

% ------------------ Optimization for Wavelet Denoising ------------------

initial_level = 3;
best_wd_level = initial_level;
min_error_wd = Inf;

for level = 1:9
    smoothed_forcing_vector_wd = wdenoise(noisy_forcing_vector, level);
    smoothed_forcing_matrix_wd = repmat(smoothed_forcing_vector_wd, 1, n_inputs);
    [y_sim_x_wd, ~] = lsim(sys_x, smoothed_forcing_matrix_wd, time_vector, xReg(L(1), 1:r-1));
    error_wd = sqrt(sum((y_sim_x_wd - y_sim_x).^2, 2));
    mean_error_wd = mean(error_wd);

    if mean_error_wd < min_error_wd
        min_error_wd = mean_error_wd;
        best_wd_level = level;
    end
end

smoothed_forcing_vector_wd_initial = wdenoise(noisy_forcing_vector, initial_level);
smoothed_forcing_vector_wd_best = wdenoise(noisy_forcing_vector, best_wd_level);

% ------------------ Plot Results ------------------

methods = {'Moving Average', 'Savitzky-Golay', 'Wavelet Denoising'};
initial_vectors = {smoothed_forcing_vector_ma_initial, smoothed_forcing_vector_sg_initial, smoothed_forcing_vector_wd_initial};
optimized_vectors = {smoothed_forcing_vector_ma_best, smoothed_forcing_vector_sg_best, smoothed_forcing_vector_wd_best};

for i = 1:3
    % Plot 1: Forcing Vector Comparison
    figure;
    plot(time_vector, forcing_vector, 'b', 'LineWidth', 1);
    hold on;
    plot(time_vector, noisy_forcing_vector, 'r', 'LineWidth', 1);
    plot(time_vector, initial_vectors{i}, 'g', 'LineWidth', 1);
    plot(time_vector, optimized_vectors{i}, 'k', 'LineWidth', 1);
    legend('Original', 'Noisy', 'Initial', 'Optimized');
    title([methods{i} ' Forcing Vector Optimization']);
    xlabel('Time');
    ylabel('Forcing Vector');
    grid on;

    % Plot 2: Error Comparison
    figure;
    plot(time_vector, error_noisy, 'r', 'LineWidth', 1);
    hold on;
    plot(time_vector, sqrt(sum((initial_vectors{i} - forcing_vector).^2, 2)), 'g', 'LineWidth', 1);
    plot(time_vector, sqrt(sum((optimized_vectors{i} - forcing_vector).^2, 2)), 'k', 'LineWidth', 1);
    legend('Noisy', 'Initial', 'Optimized');
    title([methods{i} ' Error Comparison']);
    xlabel('Time');
    ylabel('Error');
    grid on;
end

% Final Plot: Comparison of Minimum Errors
figure;
bar([mean_error_noisy, min_error_ma, min_error_sg, min_error_wd]);
set(gca, 'XTickLabel', {'Noisy', 'MA Filtered', 'SG Filtered', 'Wavelet Denoised'});
title('Comparison of Minimum Errors After Optimization');
xlabel('Method');
ylabel('Mean Error');
grid on;

% ------------------ Final Time-Series Comparison Plots ------------------

% Plot 1: Time-Series of Forcing Vectors (Original, Noisy, and Optimized)
figure;
plot(time_vector, forcing_vector, 'b', 'LineWidth', 1);
hold on;
plot(time_vector, noisy_forcing_vector, 'r', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_ma_best, 'g', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_sg_best, 'm', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_wd_best, 'k', 'LineWidth', 1);
legend('Original', 'Noisy', 'MA Optimized', 'SG Optimized', 'WD Optimized');
title('Comparison of Forcing Vectors (Original vs. Noisy vs. Optimized)');
xlabel('Time');
ylabel('Forcing Vector');
grid on;

% Plot 2: Time-Series of Errors (Noisy and Optimized Methods)
error_ma_best = sqrt(sum((smoothed_forcing_vector_ma_best - forcing_vector).^2, 2));
error_sg_best = sqrt(sum((smoothed_forcing_vector_sg_best - forcing_vector).^2, 2));
error_wd_best = sqrt(sum((smoothed_forcing_vector_wd_best - forcing_vector).^2, 2));

figure;
plot(time_vector, error_noisy, 'r', 'LineWidth', 1);
hold on;
plot(time_vector, error_ma_best, 'g', 'LineWidth', 1);
plot(time_vector, error_sg_best, 'm', 'LineWidth', 1);
plot(time_vector, error_wd_best, 'k', 'LineWidth', 1);
legend('Noisy', 'MA Optimized', 'SG Optimized', 'WD Optimized');
title('Comparison of Errors Over Time');
xlabel('Time');
ylabel('Error');
grid on;



% ------------------ End of Script ------------------
