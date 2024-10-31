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
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, sol] = ode45(lorenz, tspan, initial_conditions, options);

%Change varaible names
t2=t;
sol2=sol;
dt2=dt;

% Save the data
save('lorenzData2.mat', 'sol2', 't2', 'dt2')

% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Creating the Hankel matrix for both x and y data
m = 100;         % Size of the Hankel matrix (number of rows)
HankelMatrix_x = hankel(x(1:m), x(m:end));
HankelMatrix_y = hankel(y(1:m), y(m:end));

% Perform SVD for both x and y
[U_x, E_x, V_x] = svd(HankelMatrix_x, 'econ');
[U_y, E_y, V_y] = svd(HankelMatrix_y, 'econ');

% Derivatives computation for both x and y
r = 15;
dV_x = zeros(length(V_x)-5, r);
dV_y = zeros(length(V_y)-5, r);

for i = 3:length(V_x)-3
    for k = 1:r
        % Derivatives for V_x
        dV_x(i-2,k) = (1/(12*dt)) * (-V_x(i+2,k) + 8*V_x(i+1,k) - 8*V_x(i-1,k) + V_x(i-2,k));
        
        % Derivatives for V_y
        dV_y(i-2,k) = (1/(12*dt)) * (-V_y(i+2,k) + 8*V_y(i+1,k) - 8*V_y(i-1,k) + V_y(i-2,k));
    end
end

% Reduced variables from SVD for both x and y
xReg = V_x(3:end-3, 1:r);
yReg = V_y(3:end-3, 1:r);

% Build library of nonlinear time series for x and y
polyorder = 1;
Theta_x = poolData(xReg, r, polyorder, 0);
Theta_y = poolData(yReg, r, polyorder, 0);

% Normalize columns of Theta_x and Theta_y
for k=1:size(Theta_x,2)
    normTheta_x(k) = norm(Theta_x(:,k));
    Theta_x(:,k) = Theta_x(:,k) / normTheta_x(k);
end

for k=1:size(Theta_y,2)
    normTheta_y(k) = norm(Theta_y(:,k));
    Theta_y(:,k) = Theta_y(:,k) / normTheta_y(k);
end

% Sparse regression to obtain Xi for both x and y
lambda = 0.02;
clear Xi_x Xi_y
for k = 1:r-1
    Xi_x(:,k) = sparsifyDynamics(Theta_x, dV_x(:,k), lambda*k, 1);  
    Xi_y(:,k) = sparsifyDynamics(Theta_y, dV_y(:,k), lambda*k, 1);  
end

% Normalize Xi_x and Xi_y
for k=1:length(Xi_x)
    Xi_x(k,:) = Xi_x(k,:) / normTheta_x(k);
end

for k=1:length(Xi_y)
    Xi_y(k,:) = Xi_y(k,:) / normTheta_y(k);
end

% Extract system matrices A and B for both x and y reconstructions
A_x = Xi_x(2:r+1,1:r-1)';
B_x = A_x(:,r);
A_x = A_x(:,1:r-1);

A_y = Xi_y(2:r+1,1:r-1)';
B_y = A_y(:,r);
A_y = A_y(:,1:r-1);

%Change variable names
V_x2=V_x;
V_y2=V_y;
A_x2=A_x;
B_x2=B_x;
A_y2=A_y;
B_y2=B_y;
xReg2=xReg;
yReg2=yReg;
r2=r;
tspan2=tspan;
dt2=dt;


% Save the system data for simulation of both x and y
save('systemData2.mat', 'V_x2', 'V_y2', 'A_x2', 'B_x2', 'A_y2', 'B_y2', 'xReg2', 'yReg2', 'r2', 'tspan2', 'dt2')


