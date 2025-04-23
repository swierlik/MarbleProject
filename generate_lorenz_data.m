%generate_lorenz_data.m

% Lorenz Attractor Data Generation

% Parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Time span
dt=0.001;
tfinal=50;
tspan = [dt:dt:tfinal]; % Run for 50 time units

% Initial condition
initial_conditions = [-8, 8, 27];

% Define the Lorenz system as a function
lorenz = @(t, xyz) [sigma*(xyz(2) - xyz(1)); ...
                    xyz(1)*(rho - xyz(3)) - xyz(2); ...
                    xyz(1)*xyz(2) - beta*xyz(3)];

% Solve using ode45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, sol] = ode45(lorenz, tspan, initial_conditions, options);

% Save the data
save('Data/lorenzData.mat', 'sol', 't', 'dt')
disp('Data generation complete.')