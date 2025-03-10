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

% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Creating the Hankel matrix for x
m = 100;         % Size of the Hankel matrix (number of rows)
HankelMatrix_x = hankel(x(1:m), x(m:end));

% Perform SVD for x
[U_x, E_x, V_x] = svd(HankelMatrix_x, 'econ');

% Derivatives computation for x
r = 7;%!!!!!!!!!!!!!!!!!!!!!!
dV_x = zeros(length(V_x)-5, r);

for i = 3:length(V_x)-3
    for k = 1:r
        % Derivatives for V_x
        dV_x(i-2,k) = (1/(12*dt)) * (-V_x(i+2,k) + 8*V_x(i+1,k) - 8*V_x(i-1,k) + V_x(i-2,k));
    end
end

% Reduced variables from SVD for x
xReg = V_x(3:end-3, 1:r);

% Build library of nonlinear time series for x
polyorder = 1;
Theta_x = poolData(xReg, r, polyorder, 0);

% Normalize columns of Theta_x
for k=1:size(Theta_x,2)
    normTheta_x(k) = norm(Theta_x(:,k));
    Theta_x(:,k) = Theta_x(:,k) / normTheta_x(k);
end

% Sparse regression to obtain Xi for x
lambda = 0.0;
clear Xi_x
for k = 1:r-1
    Xi_x(:,k) = sparsifyDynamics(Theta_x, dV_x(:,k), lambda*k, 1);  
end

% Normalize Xi_x
for k=1:length(Xi_x)
    Xi_x(k,:) = Xi_x(k,:) / normTheta_x(k);
end


% Extract system matrices A and B for x reconstruction
A_x = Xi_x(2:r+1,1:r-1)';
B_x = A_x(:,r);
A_x = A_x(:,1:r-1);

% Save the system data for simulation of both x and y
save('Data/systemData.mat','V_x', 'A_x', 'B_x', 'xReg', 'tspan', 'dt', 'E_x', "U_x")

disp('Data generation complete.')


