% Load the data
load('lorenzData.mat')
load('systemData.mat')

load('lorenzData2.mat')
load('systemData2.mat')


% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
z = sol(:,3);

%UNALETERED DATA

% % Plot the original Lorenz attractor
% figure;
% subplot(1, 2, 1);
% set(gcf, 'Name', 'Lorenz Attractor');
% set(gcf, 'NumberTitle', 'off');
% plot3(x, y, z, 'LineWidth', 1.5);
% title('Lorenz Attractor');
% xlabel('X axis');
% ylabel('Y axis');
% zlabel('Z axis');
% grid on;
% view(45, 45);
% axis tight;

% % X vs Time Plot
% subplot(1, 2, 2);
% plot(t, sol(:,1), 'b', 'LineWidth', 1.5);
% title('X vs Time');
% xlabel('Time');
% ylabel('X');
% grid on;

% ----------SIMULATING THE SYSTEM--------------

% Reconstruct and simulate the system for x
L = 1:min(length(tspan), size(xReg, 1));
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Reconstruct and simulate the system for y
sys_y = ss(A_y, B_y, eye(r-1), 0*B_y);  % System matrices for y
[y_sim_y, t_sim_y] = lsim(sys_y, yReg(L, r), dt*(L-1), yReg(1, 1:r-1));

%Reconstruct and simulate the system for x with y's forcing vector
%[y_sim_x_y, t_sim_x_y] = lsim(sys_x, -1*yReg(L, r), dt*(L-1), xReg(1, 1:r-1));
[y_sim_x_y, t_sim_x_y] = lsim(sys_x, yReg(L, r), dt*(L-1), xReg(1, 1:r-1));

%Reconstruct and simulate the system for y with x's forcing vector
%[y_sim_y_x, t_sim_y_x] = lsim(sys_y, -1*xReg(L, r), dt*(L-1), yReg(1, 1:r-1));
[y_sim_y_x, t_sim_y_x] = lsim(sys_y, xReg(L, r), dt*(L-1), yReg(1, 1:r-1));

% ---------------PLOTTING X-------------------

% Plot delay-embedded attractor for x
figure;
subplot(1, 3, 1);
set(gcf, 'Name', 'Delay Embedded Attractor X and Error');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Original delay embedding of x
title('Original Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the reconstructed attractor for x
subplot(1, 3, 2);
L = 300:length(tspan)-300;
plot3(y_sim_x(L,1), y_sim_x(L,2), y_sim_x(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);

% Compute error over time
error_x = sqrt(sum((V_x(L,1:r-1) - y_sim_x(L,1:r-1)).^2, 2));
error_x2 = sqrt(sum((V_x2(L,1:r-1) - V_x(L,1:r-1)).^2, 2));

% Plot error over time
subplot(1, 3, 3);
plot(tspan(L), error_x, 'r-', 'LineWidth', 1.5);
title('Error Between Original and Reconstructed X Over Time');
xlabel('Time');
ylabel('Error');
grid on;
axis tight;

% ---------------PLOTTING Y-------------------
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

% Compute error over time
error_y = sqrt(sum((V_y(L,1:r-1) - y_sim_y(L,1:r-1)).^2, 2));
error_y2 = sqrt(sum((V_y2(L,1:r-1) - V_y(L,1:r-1)).^2, 2));



% ---------------PLOTTING X WITH Y FORCING VECTOR-------------------
% Plot delay-embedded attractor for x
figure;
subplot(1, 3, 1);
set(gcf, 'Name', 'Delay Embedded Attractor X with Y Forcing and Error');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Original delay embedding of x
title('Original Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the reconstructed attractor for x with y's forcing vector
subplot(1, 3, 2);
plot3(y_sim_x_y(L,1), y_sim_x_y(L,2), y_sim_x_y(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Attractor X with Y Forcing');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);

% Compute error over time!!!!!!!!!!!!!!!!!!!!!!!!!!
error_x_y = sqrt(sum((V_x(L,1:r-1) - y_sim_x_y(L,1:r-1)).^2, 2));
%error_x_y = sqrt(sum((V_x(L,1:3) - y_sim_x_y(L,1:3)).^2, 2));

% Plot error over time
subplot(1, 3, 3);
plot(tspan(L), error_x_y, 'r-', 'LineWidth', 1.5);
title('Error Between Original X and Reconstructed X with Y Forcing');
xlabel('Time');
ylabel('Error');
grid on;
axis tight;


% -------------COMPARING X vs XY vs X2 ERROR----------------
% plot error x vs error xy on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'Error X vs Error XY vs Error X2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L), error_x, 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L), error_x_y, 'b', 'LineWidth', 0.5)
plot(tspan(L), error_x2, 'g', 'LineWidth', 0.5)
legend("Error X", "Error XY", "Error X2")
box on;

% Define plotting range to exclude transients
L_plot = 300:length(tspan)-300;

% ---------------PLOTTING Y WITH X FORCING VECTOR-------------------
% Plot delay-embedded attractor for y
figure;
subplot(1, 3, 1);
set(gcf, 'Name', 'Delay Embedded Attractor Y with X Forcing and Error');
set(gcf, 'NumberTitle', 'off');
plot3(V_y(:,1), V_y(:,2), V_y(:,3));  % Original delay embedding of y
title('Original Delay Embedded Attractor Y');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the reconstructed attractor for y with x's forcing vector
subplot(1, 3, 2);
plot3(y_sim_y_x(L_plot,1), y_sim_y_x(L_plot,2), y_sim_y_x(L_plot,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Attractor Y with X Forcing');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);

% Compute error over time
error_y_x = sqrt(sum((V_y(L_plot,1:r-1) - y_sim_y_x(L_plot,1:r-1)).^2, 2));

% Plot error over time
subplot(1, 3, 3);
plot(tspan(L_plot), error_y_x, 'r-', 'LineWidth', 1.5);
title('Error Between Original Y and Reconstructed Y with X Forcing');
xlabel('Time');
ylabel('Error');
grid on;
axis tight;

% -------------COMPARING Y vs YX vs Y2 ERROR----------------
% plot error y vs error yx on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'Error Y vs Error YX vs Error Y2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), error_y, 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), error_y_x, 'b', 'LineWidth', 0.5)
plot(tspan(L_plot), error_y2, 'g', 'LineWidth', 0.5)
legend("Error Y", "Error YX", "Error Y2")
box on;


% -------------COMPARING X vs X2 and Y vs Y2----------------
% plot x vs x2 and y vs y2 on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'X vs X2 and Y vs Y2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_x(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_x_y(L_plot,1), 'b', 'LineWidth', 0.5)
legend("X", "X2")
box on;

figure
% Set the figure's name
set(gcf, 'Name', 'Y vs Y2 and X vs X2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_y(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_y_x(L_plot,1), 'b', 'LineWidth', 0.5)
legend("Y", "Y2")
box on;


% -------------HELPER PLOTS----------------

% % Plot the v1 of x vs v15 of y
% figure
% % Set the figure's name
% set(gcf, 'Name', 'V1x vs V15y');
% set(gcf, 'NumberTitle', 'off');
% subplot(2,1,1)
% plot(tspan(L), y_sim_x(L,1), 'b', 'LineWidth', 0.5)
% title("V1x")
% box on;

% subplot(2,1,2)
% plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
% title("V15y")
% box on;



%Plot v11 of y vs v11 of x on 1 graph
figure
% Set the figure's name
set(gcf, 'Name', 'V11y vs -V11x');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L), y_sim_x(L,11), 'b', 'LineWidth', 0.5)
%plot(tspan(L), -1*y_sim_x(L,11), 'g', 'LineWidth', 0.5)%inverse
legend("V11y", "V11x", "-V11X")
box on;


% %plot from v1 to v14 of x
% figure
% % Set the figure's name
% set(gcf, 'Name', 'V1 to V15 of x');
% set(gcf, 'NumberTitle', 'off');
% for i = 1:14
%     subplot(5,3,i)
%     plot(tspan(L), y_sim_x(L,i), 'b', 'LineWidth', 0.5)
%     title("V" + i + "x")
%     box on;
% end

% plot 