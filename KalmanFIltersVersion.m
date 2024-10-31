% Adaptive Extended Kalman Filter for Lorenz Attractor

clear all; close all; clc;

%% Parameters of the Lorenz system
sigma = 10;
beta = 8/3;
rho = 28;

% Time parameters
dt = 0.001;
T = 50;
t = 0:dt:T;
N = length(t);

%% Initialize true state
x_true = zeros(3,N);
x_true(:,1) = [-8; 8; 27];  % Initial condition

%% Simulate the Lorenz system
for k = 1:N-1
    x = x_true(:,k);
    x_true(:,k+1) = x + lorenz_system(x, sigma, beta, rho)*dt;
end

%% Add measurement noise
% Assume we can measure all states with noise
measurement_noise_std = [2; 2; 2];
z = x_true + measurement_noise_std .* randn(3,N);

%% Adaptive Extended Kalman Filter Initialization
% Initial estimate
x_hat = zeros(3,N);
x_hat(:,1) = [0; 0; 0];  % Initial guess

% Initial estimation error covariance
P = eye(3);

% Process noise covariance (tuned adaptively)
Q = eye(3);

% Measurement noise covariance
R = diag(measurement_noise_std.^2);

% Identity matrix
I = eye(3);

%% Adaptive EKF Implementation
for k = 1:N-1
    % Prediction Step
    x_pred = x_hat(:,k) + lorenz_system(x_hat(:,k), sigma, beta, rho)*dt;
    
    % Jacobian of the system
    F = jacobian_lorenz(x_hat(:,k), sigma, beta, rho)*dt + I;
    
    % Predicted estimation error covariance
    P_pred = F * P * F' + Q;
    
    % Measurement Matrix (Identity, since we measure all states)
    H = I;
    
    % Innovation
    y = z(:,k+1) - x_pred;
    
    % Innovation covariance
    S = H * P_pred * H' + R;
    
    % Kalman Gain
    K = P_pred * H' / S;
    
    % Update Step
    x_hat(:,k+1) = x_pred + K * y;
    
    % Update estimation error covariance
    P = (I - K * H) * P_pred;
    
    % Adaptive Tuning of Q
    Q = adapt_process_noise(Q, y, S, K);
end

%% Plot Results
figure;
plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'b', 'LineWidth', 1.5); hold on;
plot3(x_hat(1,:), x_hat(2,:), x_hat(3,:), 'r--', 'LineWidth', 1.5);
legend('True State', 'Estimated State');
xlabel('x'); ylabel('y'); zlabel('z');
title('Lorenz Attractor State Estimation using Adaptive EKF');
grid on;


%%Plot error
error = sqrt(sum((x_true - x_hat).^2, 1));

figure;
plot(t, error, 'r-', 'LineWidth', 1.5);
title('Error Between True and Estimated State');
xlabel('Time');
ylabel('Error');
grid on;

%% Supporting Functions

function dx = lorenz_system(x, sigma, beta, rho)
    % Lorenz system differential equations
    dx = zeros(3,1);
    dx(1) = sigma * (x(2) - x(1));
    dx(2) = x(1) * (rho - x(3)) - x(2);
    dx(3) = x(1) * x(2) - beta * x(3);
end

function F = jacobian_lorenz(x, sigma, beta, rho)
    % Jacobian matrix of the Lorenz system
    F = zeros(3,3);
    F(1,1) = -sigma;
    F(1,2) = sigma;
    F(2,1) = rho - x(3);
    F(2,2) = -1;
    F(2,3) = -x(1);
    F(3,1) = x(2);
    F(3,2) = x(1);
    F(3,3) = -beta;
end

function Q_new = adapt_process_noise(Q_prev, y, S, K)
    % Adaptively update the process noise covariance Q
    % This is a simple heuristic adaptation
    alpha = 0.01;  % Learning rate
    e = y - K * y;  % Estimation error
    Q_new = Q_prev + alpha * (K * e * e' * K');
end
