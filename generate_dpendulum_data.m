%generate_dpendulum_data.m
% This script generates data for the double pendulum system and saves it to a .mat file.

%Parameters
g = 9.81; % Acceleration due to gravity
l1 = 1; % Length of the first pendulum
l2 = 1; % Length of the second pendulum
m1 = 1; % Mass of the first pendulum
m2 = 1; % Mass of the second pendulum

% Time span
dt = 0.001; % Time step
tfinal = 50; % Final time
tspan = 0:dt:tfinal; % Time vector

% Initial conditions
theta1_0 = pi/4; % Initial angle of the first pendulum
theta2_0 = pi/4; % Initial angle of the second pendulum
theta1_dot_0 = 0; % Initial angular velocity of the first pendulum
theta2_dot_0 = 0; % Initial angular velocity of the second pendulum

initial_conditions = [theta1_0; theta1_dot_0; theta2_0; theta2_dot_0]; % Initial state vector

dpendulum = @(t, y) [y(2); ... % Angular velocity of the first pendulum
                    (-g*(2*m1 + m2)*sin(y(1)) - m2*g*sin(y(1) - 2*y(3)) - 2*sin(y(1) - y(3))*m2*(y(4)^2*l2 + y(2)^2*l1*cos(y(1) - y(3)))) / (l1*(2*m1 + m2 - m2*cos(2*y(1) - 2*y(3)))); ... % Angular acceleration of the first pendulum
                    y(4); ... % Angular velocity of the second pendulum
                    (2*sin(y(1) - y(3))*(y(2)^2*l1*(m1 + m2) + g*(m1 + m2)*cos(y(1)) + y(4)^2*l2*m2*cos(y(1) - y(3)))) / (l2*(m1 + m2 - m2*cos(2*y(1) - 2*y(3)))); ... % Angular acceleration of the second pendulum
                    ];

%Solve using ode45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % Set ODE solver options
[t, sol] = ode45(dpendulum, tspan, initial_conditions, options); % Solve the ODE

save('Data/dpendulumData.mat', 'sol', 't', 'dt') % Save the data to a .mat file
disp('Data generation complete.') % Display completion message