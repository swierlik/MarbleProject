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


% Reconstruct and simulate the system for x
L = 1:min(50000, size(xReg, 1));
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Reconstruct and simulate the system for y
sys_y = ss(A_y, B_y, eye(r-1), 0*B_y);  % System matrices for y
[y_sim_y, t_sim_y] = lsim(sys_y, yReg(L, r), dt*(L-1), yReg(1, 1:r-1));

%Reconstruct and simulate the system for x with y's forcing vector
sys_x_y = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x_y, t_sim_x_y] = lsim(sys_x_y, yReg(L, r), dt*(L-1), xReg(1, 1:r-1));

save("reconstructedData.mat", "y_sim_x", "y_sim_y", "y_sim_x_y", "t_sim_x", "t_sim_y", "t_sim_x_y")

%load("reconstructedData.mat")


%PLOTTING X
% Plot delay-embedded attractor for x
figure;
subplot(1, 2, 1);
set(gcf, 'Name', 'Delay Embedded Attractor');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Using V_x for the delay embedding of x
title('Delay Embedded Attractor X');
view(-15, 65);

% Plot the reconstructed attractor for x
subplot(1, 2, 2);
set(gcf, 'Name', 'Reconstructed Delay Embedded Attractor X');
set(gcf, 'NumberTitle', 'off');
L = 300:length(tspan)-300;
plot3(y_sim_x(L,1), y_sim_x(L,2), y_sim_x(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor X');
axis tight;
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

%PLOTTING Y
% Plot delay-embedded attractor for y
figure;
subplot(1, 2, 1);
set(gcf, 'Name', 'Delay Embedded Attractor');
set(gcf, 'NumberTitle', 'off');
plot3(V_y(:,1), V_y(:,2), V_y(:,3));  % Using V_y for the delay embedding of y
title('Delay Embedded Attractor Y');
view(-15, 65);

% Plot the reconstructed attractor for y
subplot(1, 2, 2);
set(gcf, 'Name', 'Reconstructed Delay Embedded Attractor Y');
set(gcf, 'NumberTitle', 'off');
L = 300:length(tspan)-300;
plot3(y_sim_y(L,1), y_sim_y(L,2), y_sim_y(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor Y');
axis tight;
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

%PLOTTING XY
% Plot delay-embedded attractor for x
figure;
subplot(1, 2, 1);
set(gcf, 'Name', 'Delay Embedded Attractor');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Using V_x for the delay embedding of x
title('Delay Embedded Attractor X');
view(-15, 65);

% Plot the reconstructed attractor for x with y's forcing vector
subplot(1, 2, 2);
set(gcf, 'Name', 'Reconstructed Delay Embedded Attractor X with Y forcing vector');
set(gcf, 'NumberTitle', 'off');
L = 300:length(tspan)-300;
plot3(y_sim_x_y(L,1), y_sim_x_y(L,2), y_sim_x_y(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor X with Y forcing vector');
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
title("V1x")
box on;

subplot(2,1,2)
plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
title("V15y")
box on;

%Plot v11 of y vs v11 of x on 1 graph
figure
% Set the figure's name
set(gcf, 'Name', 'V11y vs V11x');
% Optional: Make the name visible in the figure window title bar
set(gcf, 'NumberTitle', 'off');
plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L), y_sim_x(L,11), 'b', 'LineWidth', 0.5)
legend("V11y", "V11x")
box on;



%plot v1 to v15 of x
figure
% Set the figure's name
set(gcf, 'Name', 'V1 to V15 of x');
% Optional: Make the name visible in the figure window title bar
set(gcf, 'NumberTitle', 'off');
for i = 1:15
    subplot(5,3,i)
    plot(tspan(L), y_sim_x(L,i), 'b', 'LineWidth', 0.5)
    title("V" + i + "x")
    box on;
end
