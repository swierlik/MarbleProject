% Lorenz Attractor

% Parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Time span
dt=0.001;
tspan = [dt:dt:200]; % Run for 250 time units

% Initial condition
initial_conditions = [-8, 8, 27];

% Define the Lorenz system as a function
lorenz = @(t, xyz) [sigma*(xyz(2) - xyz(1)); ...
                    xyz(1)*(rho - xyz(3)) - xyz(2); ...
                    xyz(1)*xyz(2) - beta*xyz(3)];

% Solve using ode45
[t, sol] = ode45(lorenz, tspan, initial_conditions);

% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Plot the Lorenz attractor
figure;
plot3(x, y, z, 'LineWidth', 1.5);
title('Lorenz Attractor');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
grid on;
view(45,45);  % Set the view to 3D
axis tight;


% Creating the Hankel matrix
m = 100;         % Size of the Hankel matrix (number of rows)
HankelMatrix = hankel(x(1:m),x(m:end));

[U,E,V]=svd(HankelMatrix,'econ');

figure;
plot3(V(:,1),V(:,2),V(:,3))

% Next step aka derivatives
r=15;

dV = zeros(length(V)-5, r);
for i = 3:length(V)-3
    for k = 1:r
        dV(i-2,k) = (1/(12*dt)) * (-V(i+2,k) + 8*V(i+1,k) - 8*V(i-1,k) + V(i-2,k));
    end
end

xReg = V(3:end-3, 1:r);
dx = dV;

% Build library of nonlinear time series
polyorder = 1;
Theta = poolData(xReg,r,1,0);
% normalize columns of Theta (required in new time-delay coords)
for k=1:size(Theta,2)
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
m = size(Theta,2);
% compute Sparse regression: sequential least squares
% requires different lambda parameters for each column
lambda=0;
clear Xi
for k=1:r-1
    Xi(:,k) = sparsifyDynamics(Theta,dx(:,k),lambda*k,1);  % lambda = 0 gives better results 
end
Theta = poolData(xReg,r,1,0);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end
A = Xi(2:r+1,1:r-1)';
B = A(:,r);
A = A(:,1:r-1);
% Simulate System
L = 1:50000;
sys = ss(A, B, eye(r-1), 0*B);
[y_sim, t_sim] = lsim(sys, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

figure
L = 300:50000;
plot3(y_sim(L,1),y_sim(L,2),y_sim(L,3),'Color',[0 0 .5],'LineWidth',1.5)
axis tight
xlabel('v_1'), ylabel('v_2'), zlabel('v_3')