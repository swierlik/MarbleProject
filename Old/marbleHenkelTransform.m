% Lorenz Attractor

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


%Save the data
%save('lorenzData.mat','sol','t','dt')

% Load the data
%load('lorenzData.mat')

% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);


% Creating the Hankel matrix
m = 100;         % Size of the Hankel matrix (number of rows)
HankelMatrix = hankel(x(1:m),x(m:end));

%add some noise to hankel
HankelMatrix = HankelMatrix + 0.01*randn(size(HankelMatrix));


[U,E,V]=svd(HankelMatrix,'econ');

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

%Sparse
lambda=0.02;
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
disp(A)
% Simulate System
L = 1:min(50000, size(xReg, 1));
sys = ss(A, B, eye(r-1), 0*B);
[y_sim, t_sim] = lsim(sys, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));




% % Plot the Lorenz attractor
% figure;
% subplot(1, 2, 1);
% % Set the figure's name
% set(gcf, 'Name', 'Lorenz Attractor');
% % Optional: Make the name visible in the figure window title bar
% set(gcf, 'NumberTitle', 'off');
% plot3(x, y, z, 'LineWidth', 1.5);
% title('Lorenz Attractor');
% xlabel('X axis');
% ylabel('Y axis');
% zlabel('Z axis');
% grid on;
% view(45,45);  % Set the view to 3D
% axis tight;

% % X vs Time Plot
% subplot(1, 2, 2);
% plot(t, sol(:,1), 'b', 'LineWidth', 1.5);
% title('X vs Time');
% xlabel('Time');
% ylabel('X');
% grid on;

% figure
% % Set the figure's name
% set(gcf, 'Name', 'Delay Embedded Attractor');
% % Optional: Make the name visible in the figure window title bar
% set(gcf, 'NumberTitle', 'off');
% plot3(V(:,1),V(:,2),V(:,3))
% view(-15,65);  % Set the view to 3D

figure
% Set the figure's name
set(gcf, 'Name', 'Reconstructed Attractor with noise');
% Optional: Make the name visible in the figure window title bar
set(gcf, 'NumberTitle', 'off');
L = 300:length(tspan)-300;
plot3(y_sim(L,1),y_sim(L,2),y_sim(L,3),'Color',[0 0 .5],'LineWidth',1.5)
axis tight
xlabel('v_1'), ylabel('v_2'), zlabel('v_3')
view(-15,65);  % Set the view to 3D


% figure
% % Set the figure's name
% set(gcf, 'Name', 'V1 vs V15');
% % Optional: Make the name visible in the figure window title bar
% set(gcf, 'NumberTitle', 'off');
% subplot(2,1,1)
% plot(tspan(L), y_sim(L,1), 'b', 'LineWidth', 0.5)
% box on;

% subplot(2,1,2)
% plot(tspan(L), xReg(L,11), 'r', 'LineWidth', 0.5)
% box on;

% figure
% % Set the figure's name
% set(gcf, 'Name', 'Original vs Reconstructed');
% % Optional: Make the name visible in the figure window title bar
% set(gcf, 'NumberTitle', 'off');
% plot(tspan(L), sol(L,1), 'r', 'LineWidth', 0.5)
% hold on
% plot(tspan(L), y_sim(L,1), 'b', 'LineWidth', 0.5)


% %%Recurrence plot
% figure
% % Set the figure's name
% set(gcf, 'Name', 'Recurrence Plot');
% % Optional: Make the name visible in the figure window title bar
% set(gcf, 'NumberTitle', 'off');
% % Recurrence plot
% %imagesc(abs(xcorr(x,x)>0.1));
% % Compute the pairwise distances between points
% n = length(x);
% RP = zeros(n);  % Initialize the recurrence matrix (this will store distances)

% for i = 1:n
%     for j = 1:n
%         RP(i, j) = abs(x(i) - x(j));  % Store the distance between points
%     end
% end
% imagesc(RP);
% title('Recurrence Plot');
% xlabel('Time');
% ylabel('Time');
% colormap('gray');
% axis square;
