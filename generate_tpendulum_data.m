%generate_tpendulum_data.m
% This script generates data for the triple pendulum system and saves it to a .mat file.

%Parameters
g = 9.81; % Acceleration due to gravity
l1 = 1; % Length of the first pendulum
l2 = 1; % Length of the second pendulum
l3 = 1; % Length of the third pendulum
m1 = 1; % Mass of the first pendulum
m2 = 1; % Mass of the second pendulum
m3 = 1; % Mass of the third pendulum

% Time span
dt = 0.001; % Time step
tfinal = 50; % Final time
tspan = 0:dt:tfinal; % Time vector

%initial conditions
theta1_0 = pi/4; % Initial angle of the first pendulum
theta2_0 = pi/4; % Initial angle of the second pendulum
theta3_0 = pi/4; % Initial angle of the third pendulum
theta1_dot_0 = 0; % Initial angular velocity of the first pendulum
theta2_dot_0 = 0; % Initial angular velocity of the second pendulum
theta3_dot_0 = 0; % Initial angular velocity of the third pendulum

initial_conditions = [theta1_0; theta1_dot_0; theta2_0; theta2_dot_0; theta3_0; theta3_dot_0]; % Initial state vector

tpendulum = @(t, y) [y(2); ... % Angular velocity of the first pendulum
                     (-g*(2*m1 + m2 + m3)*sin(y(1)) - m2*g*sin(y(1) - 2*y(3)) - m3*g*sin(y(1) - 2*y(5)) - 2*sin(y(1) - y(3))*m2*(y(4)^2*l2 + y(2)^2*l1*cos(y(1) - y(3))) - 2*sin(y(1) - y(5))*m3*(y(6)^2*l3 + y(2)^2*l1*cos(y(1) - y(5)))) / (l1*(2*m1 + m2 + m3 - m2*cos(2*y(1) - 2*y(3)) - m3*cos(2*y(1) - 2*y(5)))); ... % Angular acceleration of the first pendulum
                     y(4); ... % Angular velocity of the second pendulum
                     (g*(m1 + m2 + m3)*cos(y(1)) + y(4)^2*l2*m2*cos(y(1) - y(3)) + y(6)^2*l3*m3*cos(y(1) - y(5))) / (l2*(m1 + m2 + m3 - m2*cos(2*y(1) - 2*y(3)) - m3*cos(2*y(1) - 2*y(5)))); ... % Angular acceleration of the second pendulum
                     y(6); ... % Angular velocity of the third pendulum
                     (g*(m1 + m2 + m3)*cos(y(1)) + y(4)^2*l2*m2*cos(y(1) - y(5)) + y(6)^2*l3*m3*cos(y(1) - y(5))) / (l3*(m1 + m2 + m3 - m3*cos(2*y(1) - 2*y(5)))); ... % Angular acceleration of the third pendulum
                     ];

%Solve using ode45
%options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8); % Set ODE solver options

[t, sol] = ode45(tpendulum, tspan, initial_conditions); % Solve the ODE

save('Data/tpendulumData.mat', 'sol', 't', 'dt') % Save the data to a .mat file
disp('Data generation complete.') % Display completion message