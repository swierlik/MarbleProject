% noise_reduction.m

% Load the data
load('lorenzData.mat') % Contains 'sol', 't', 'dt'
load('systemData.mat') % Contains 'V_x', 'V_y', 'A_x', 'B_x', 'A_y', 'B_y', 'xReg', 'yReg', 'r', 'tspan', 'dt')

% Extract x, y, z from sol
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Define L (exclude initial and final transients)
L = 300:length(xReg)-300;

% Forcing vector (original)
forcing_vector = xReg(L, r);

% Determine the number of input channels
n_inputs = size(B_x, 2);

% Define the time vector
time_vector = dt * (L - 1);

% Simulate original system (without noise) for comparison
sys_x = ss(A_x, B_x, eye(r-1), zeros(r-1, n_inputs));
forcing_matrix = repmat(forcing_vector, 1, n_inputs);
[y_sim_x, t_sim_x] = lsim(sys_x, forcing_matrix, time_vector, xReg(L(1), 1:r-1));

% ------------------ Add Noise to the Forcing Vector ------------------
% Set noise level (adjust as needed)
noise_level = 0.1 * std(forcing_vector);

% Add Gaussian noise
noisy_forcing_vector = forcing_vector + noise_level * randn(size(forcing_vector));
noisy_forcing_matrix = repmat(noisy_forcing_vector, 1, n_inputs);

% Reconstruct and simulate the system with noisy forcing vector
[y_sim_x_noisy, t_sim_x_noisy] = lsim(sys_x, noisy_forcing_matrix, time_vector, xReg(L(1), 1:r-1));
error_noisy = sqrt(sum((y_sim_x_noisy - y_sim_x).^2, 2));

% ------------------ Apply Noise Reduction Methods ------------------

% Method 1: Moving Average Filter
window_size = 5;
smoothed_forcing_vector_ma = movmean(noisy_forcing_vector, window_size);
smoothed_forcing_matrix_ma = repmat(smoothed_forcing_vector_ma, 1, n_inputs);
[y_sim_x_ma, t_sim_x_ma] = lsim(sys_x, smoothed_forcing_matrix_ma, time_vector, xReg(L(1), 1:r-1));
error_ma = sqrt(sum((y_sim_x_ma - y_sim_x).^2, 2));

% Method 2: Savitzky-Golay Filter
order = 3;
framelen = 7;
smoothed_forcing_vector_sg = sgolayfilt(noisy_forcing_vector, order, framelen);
smoothed_forcing_matrix_sg = repmat(smoothed_forcing_vector_sg, 1, n_inputs);
[y_sim_x_sg, t_sim_x_sg] = lsim(sys_x, smoothed_forcing_matrix_sg, time_vector, xReg(L(1), 1:r-1));
error_sg = sqrt(sum((y_sim_x_sg - y_sim_x).^2, 2));

% Method 3: Wavelet Denoising
smoothed_forcing_vector_wd = wdenoise(noisy_forcing_vector);
smoothed_forcing_matrix_wd = repmat(smoothed_forcing_vector_wd, 1, n_inputs);
[y_sim_x_wd, t_sim_x_wd] = lsim(sys_x, smoothed_forcing_matrix_wd, time_vector, xReg(L(1), 1:r-1));
error_wd = sqrt(sum((y_sim_x_wd - y_sim_x).^2, 2));

% ------------------ Plotting Original Figures ------------------

methods = {'Moving Average', 'Savitzky-Golay', 'Wavelet Denoising'};
y_sims = {y_sim_x_ma, y_sim_x_sg, y_sim_x_wd};
errors = {error_ma, error_sg, error_wd};

for i = 1:3
    figure;
    subplot(1,2,1);
    plot3(y_sim_x(:,1), y_sim_x(:,2), y_sim_x(:,3), 'b', 'LineWidth', 1);
    hold on;
    plot3(y_sims{i}(:,1), y_sims{i}(:,2), y_sims{i}(:,3), 'r', 'LineWidth', 1);
    title(['Reconstructed Attractor with ' methods{i}]);
    legend('Original', methods{i});
    xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
    grid on;
    view(-15, 65);

    subplot(1,2,2);
    plot(time_vector, errors{i}, 'k', 'LineWidth', 1);
    title(['Error with ' methods{i}]);
    xlabel('Time');
    ylabel('Error');
    grid on;
end

% ------------------ New Combined Forcing Vector Plot ------------------

figure;
plot(time_vector, forcing_vector, 'b', 'LineWidth', 1);
hold on;
plot(time_vector, noisy_forcing_vector, 'r', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_ma, 'g', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_sg, 'm', 'LineWidth', 1);
plot(time_vector, smoothed_forcing_vector_wd, 'k', 'LineWidth', 1);
legend('Original', 'Noisy', 'MA Filtered', 'SG Filtered', 'Wavelet Denoised');
title('Comparison of Forcing Vectors');
xlabel('Time');
ylabel('Forcing Vector');
grid on;

% ------------------ New Combined Reconstructed Attractor Plot (v_1 Component) ------------------

figure;
plot(time_vector, y_sim_x(:,1), 'b', 'LineWidth', 1);
hold on;
plot(time_vector, y_sim_x_noisy(:,1), 'r', 'LineWidth', 1);
plot(time_vector, y_sim_x_ma(:,1), 'g', 'LineWidth', 1);
plot(time_vector, y_sim_x_sg(:,1), 'm', 'LineWidth', 1);
plot(time_vector, y_sim_x_wd(:,1), 'k', 'LineWidth', 1);
legend('Original', 'Noisy', 'MA Filtered', 'SG Filtered', 'Wavelet Denoised');
title('Comparison of Reconstructed Attractors (v_1 Component)');
xlabel('Time');
ylabel('v_1 Component');
grid on;

% ------------------ New Combined Error Plot ------------------

figure;
plot(time_vector, error_noisy, 'r', 'LineWidth', 1);
hold on;
plot(time_vector, error_ma, 'g', 'LineWidth', 1);
plot(time_vector, error_sg, 'm', 'LineWidth', 1);
plot(time_vector, error_wd, 'k', 'LineWidth', 1);
legend('Noisy', 'MA Filtered', 'SG Filtered', 'Wavelet Denoised');
title('Comparison of Errors for Different Noise Reduction Methods');
xlabel('Time');
ylabel('Error');
grid on;

% ------------------ End of Script ------------------
