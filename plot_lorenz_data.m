% Load the data
load('lorenzData.mat')
load('systemData.mat')

% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

% Plot the original Lorenz attractor
figure;
subplot(1, 2, 1);
set(gcf, 'Name', 'Lorenz Attractor');
set(gcf, 'NumberTitle', 'off');
plot3(x, y, z, 'LineWidth', 1.5);
title('Lorenz Attractor');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
grid on;
view(45, 45);
axis tight;

% X vs Time Plot
subplot(1, 2, 2);
plot(t, sol(:,1), 'b', 'LineWidth', 1.5);
title('X vs Time');
xlabel('Time');
ylabel('X');
grid on;

% Plot delay-embedded attractor for x
figure;
set(gcf, 'Name', 'Delay Embedded Attractor');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Using V_x for the delay embedding of x
view(-15, 65);

% Reconstruct and simulate the system for x
L = 1:min(50000, size(xReg, 1));
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Plot the reconstructed attractor for x
figure;
set(gcf, 'Name', 'Reconstructed Attractor');
set(gcf, 'NumberTitle', 'off');
L = 300:length(tspan)-300;
plot3(y_sim_x(L,1), y_sim_x(L,2), y_sim_x(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
axis tight;
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the v1 of x vs v15 of y
figure
% Set the figure's name
set(gcf, 'Name', 'V1x vs V15y');
% Optional: Make the name visible in the figure window title bar
set(gcf, 'NumberTitle', 'off');
subplot(2,1,1)
plot(tspan(L), y_sim_x(L,1), 'b', 'LineWidth', 0.5)
box on;

subplot(2,1,2)
plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
box on;
