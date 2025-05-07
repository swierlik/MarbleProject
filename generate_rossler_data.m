%generate_rossler_data.m

% Rössler Attractor Data Generation
%Parameters
sigma = 0.2;
a = 0.2;
b = 0.2;
c = 5.7;

%Time span
dt=0.001;
tfinal=50;
tspan = [dt:dt:tfinal]; % Run for 50 time units

%Initial condition
initial_conditions = [1, 1, 1];


%Define the Rössler system as a function
rossler = @(t, xyz) [-xyz(2) - xyz(3); ... % dx/dt
                     xyz(1) + a*xyz(2); ... % dy/dt
                     b + xyz(3)*(xyz(1) - c)]; % dz/dt

%Solve using ode45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, sol] = ode45(rossler, tspan, initial_conditions, options);

% Save the data
save('Data/rosslerData.mat', 'sol', 't', 'dt')
disp('Data generation complete.')
