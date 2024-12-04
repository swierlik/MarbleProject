% havok_reconstruction.m

% Load the data
load('Data/lorenzData.mat') % Original deterministic data
load('Data/lorenzDataStochastic.mat') % Stochastic noisy data

% Extract x, y, z components
x_original = sol(:,1); % Deterministic data
y_original = sol(:,2);
z_original = sol(:,3);

x_noisy = sol2(:,1); % Stochastic noisy data
y_noisy = sol2(:,2);
z_noisy = sol2(:,3);

% Time vector
dt = 0.001; % Ensure it matches your original data
time_vector = dt * (1:length(x_original));

% ------------------ Noise Reduction on Stochastic Data ------------------

% Universal parameters (from optimization results)
window_size = 11; % Moving Average
order = 2; % Savitzky-Golay Order
framelen = 5; % Savitzky-Golay Frame Length
wd_level = 5; % Wavelet Denoising Level

% Apply noise reduction to x, y, z
x_smoothed = movmean(x_noisy, window_size); % Moving Average (example)
y_smoothed = sgolayfilt(y_noisy, order, framelen); % Savitzky-Golay
z_smoothed = wdenoise(z_noisy, wd_level); % Wavelet Denoising

% ------------------ Generate HAVOK System Data ------------------

% Combine smoothed variables
sol_smoothed = [x_smoothed, y_smoothed, z_smoothed];

% Creating the Hankel matrix for the x data
m = 100; % Size of the Hankel matrix (number of rows)
HankelMatrix_x = hankel(x_smoothed(1:m), x_smoothed(m:end));

% Perform Singular Value Decomposition (SVD)
[U_x, E_x, V_x] = svd(HankelMatrix_x, 'econ');

% Derivative computation for HAVOK
r = 15; % Truncation rank for HAVOK
dV_x = zeros(length(V_x)-5, r);

for i = 3:length(V_x)-3
    for k = 1:r
        % Derivative for V_x
        dV_x(i-2, k) = (1/(12*dt)) * (-V_x(i+2, k) + 8*V_x(i+1, k) - 8*V_x(i-1, k) + V_x(i-2, k));
    end
end

% Reduced variables from SVD
xReg = V_x(3:end-3, 1:r);

% Build library of nonlinear time series for x
polyorder = 1; % Only linear terms
Theta_x = poolData(xReg, r, polyorder, 0);

% Normalize columns of Theta_x
for k = 1:size(Theta_x, 2)
    normTheta_x(k) = norm(Theta_x(:, k));
    Theta_x(:, k) = Theta_x(:, k) / normTheta_x(k);
end

% Sparse regression to obtain Xi for x
lambda = 0.0; % Regularization parameter
Xi_x = sparsifyDynamics(Theta_x, dV_x(:, 1:r), lambda, 1);

% Extract system matrices A and B
A_x = Xi_x(2:r+1, 1:r-1)';
B_x = A_x(:, r);
A_x = A_x(:, 1:r-1);

% ------------------ Simulate HAVOK System ------------------

L = 1:min(length(time_vector), size(xReg, 1));
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x); % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% ------------------ Plot Results ------------------

% Original vs Reconstructed Attractor
figure;
plot3(x_original, y_original, z_original, 'b', 'LineWidth', 1); % Original
hold on;
plot3(y_sim_x(:,1), y_sim_x(:,2), y_sim_x(:,3), 'r', 'LineWidth', 1); % Reconstructed
title('Original vs Reconstructed Attractor');
legend('Original', 'Reconstructed');
xlabel('x'), ylabel('y'), zlabel('z');
grid on;
view(-15, 65);

% Comparison of x-component
figure;
plot(time_vector, x_original, 'b', 'LineWidth', 1); % Original
hold on;
plot(time_vector(1:length(y_sim_x(:,1))), y_sim_x(:,1), 'r', 'LineWidth', 1); % Reconstructed
title('Comparison of x-component');
legend('Original', 'Reconstructed');
xlabel('Time');
ylabel('x');
grid on;

% Reconstruction Error
error_reconstruction = sqrt(sum((y_sim_x - [x_original(1:length(y_sim_x)), ...
                                            y_original(1:length(y_sim_x)), ...
                                            z_original(1:length(y_sim_x))]).^2, 2));

figure;
plot(time_vector(1:length(error_reconstruction)), error_reconstruction, 'k', 'LineWidth', 1);
title('Reconstruction Error Over Time');
xlabel('Time');
ylabel('Error');
grid on;

% ------------------ End of Script ------------------
