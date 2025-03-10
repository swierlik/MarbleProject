% generate_noisy_data.m

% Load the original deterministic Lorenz data
load('Data/lorenzData.mat', 'sol', 't', 'dt');

% Set Gaussian noise level
eta = 0.00000001;  % Standard deviation of Gaussian noise

% Generate measurement noise (Gaussian with mean 0 and std dev eta)
measurement_noise_x = eta * randn(size(sol,1), 1);
measurement_noise_y = eta * randn(size(sol,1), 1);
measurement_noise_z = eta * randn(size(sol,1), 1);

% Create a noisy version of the solution (simulating measurement noise)
sol2 = sol;
sol2(:, 1) = sol2(:, 1) + measurement_noise_x;
sol2(:, 2) = sol2(:, 2) + measurement_noise_y;
sol2(:, 3) = sol2(:, 3) + measurement_noise_z;

% Save the noisy dataset
save('Data/lorenzDataStochastic.mat', 'sol2', 't', 'dt');

% Extract x, y, z from noisy data
x = sol2(:, 1);
y = sol2(:, 2);
z = sol2(:, 3);

% ------------------ Build System Data ------------------

% Creating the Hankel matrix for x data
m = 100;  % Size of the Hankel matrix (number of rows)
HankelMatrix_x = hankel(x(1:m), x(m:end));

% Perform SVD for both x and y
[U_x, E_x, V_x] = svd(HankelMatrix_x, 'econ');

% Derivatives computation for both x and y
r = 15;%!!!!!!!!!!!!!!!!!!!!!!!!
dV_x = zeros(length(V_x)-5, r);


for i = 3:length(V_x)-3
    for k = 1:r
        % Derivatives for V_x
        dV_x(i-2, k) = (1/(12*dt)) * (-V_x(i+2, k) + 8*V_x(i+1, k) - 8*V_x(i-1, k) + V_x(i-2, k));
    end
end

% Reduced variables from SVD x
xReg = V_x(3:end-3, 1:r);


% Build library of nonlinear time series for x
polyorder = 1;
Theta_x = poolData(xReg, r, polyorder, 0);

% Normalize columns of Theta_x and Theta_y
for k = 1:size(Theta_x, 2)
    normTheta_x(k) = norm(Theta_x(:, k));
    Theta_x(:, k) = Theta_x(:, k) / normTheta_x(k);
end
% Sparse regression to obtain Xi for x
lambda = 0.0;
clear Xi_x
for k = 1:r-1
    Xi_x(:, k) = sparsifyDynamics(Theta_x, dV_x(:, k), lambda * k, 1);   
end

% Normalize Xi_x and Xi_y
for k = 1:length(Xi_x)
    Xi_x(k, :) = Xi_x(k, :) / normTheta_x(k);
end

% Extract system matrices A and B for both x and y reconstructions
A_x = Xi_x(2:r+1, 1:r-1)';
B_x = A_x(:, r);
A_x = A_x(:, 1:r-1);


% Change variable names for saving
V_x2 = V_x;
A_x2 = A_x;
B_x2 = B_x;
xReg2 = xReg;
r2 = r;
tspan2 = t;
dt2 = dt;
U_x2 = U_x;
E_x2 = E_x;

% Save the system data for simulation of both x and y
save('Data/systemDataStochastic.mat', 'V_x2', 'A_x2', 'B_x2', 'xReg2', 'r2', 'tspan2', 'dt2', 'U_x2', 'E_x2');

% Display completion message with noise level
disp('Stochastic data generation complete: Gaussian noise added to Lorenz data.');
disp(['Gaussian noise level (standard deviation): ', num2str(eta)]);

