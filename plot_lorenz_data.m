%plot_lorenz_data.m

% Load the data

%Lorenz data
load('Data/lorenzData.mat')
load('Data/systemData.mat')

load('Data/lorenzDataStochastic.mat')
load('Data/systemDataStochastic.mat')

% %Osillator data
% load('oscillatorData.mat')
% load('systemData.mat')

% load('oscillatorData2.mat')
% load('systemData2.mat')


% Extract x, y, z
x = sol(:,1);
y = sol(:,2);
%z = sol(:,3);

%UNALETERED DATA

% Plot the original Lorenz attractor
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
L2 = 1:min(length(tspan), size(xReg2, 1));
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Reconstruct and simulate the system for x with 2nd data
sys_x2 = ss(A_x2, B_x2, eye(r2-1), 0*B_x2);  % System matrices for x
[y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L2, r2), dt2*(L2-1), xReg2(1, 1:r2-1));

%Reconstruct and simulate the system for y with 2nd data
sys_y2 = ss(A_y2, B_y2, eye(r2-1), 0*B_y2);  % System matrices for y
[y_sim_y2, t_sim_y2] = lsim(sys_y2, yReg2(L2, r2), dt2*(L2-1), yReg2(1, 1:r2-1));

% Reconstruct and simulate the system for y
sys_y = ss(A_y, B_y, eye(r-1), 0*B_y);  % System matrices for y
[y_sim_y, t_sim_y] = lsim(sys_y, yReg(L, r), dt*(L-1), yReg(1, 1:r-1));

%Reconstruct and simulate the system for x with y's forcing vector!!!!!!!!XY=0
[y_sim_x_y, t_sim_x_y] = lsim(sys_x, yReg(L, r), dt*(L-1), xReg(1, 1:r-1));

%Reconstruct and simulate the system for y with x's forcing vector
%[y_sim_y_x, t_sim_y_x] = lsim(sys_y, -1*xReg(L, r), dt*(L-1), yReg(1, 1:r-1));
[y_sim_y_x, t_sim_y_x] = lsim(sys_y, xReg(L, r), dt*(L-1), yReg(1, 1:r-1));

%Recostuct and simulate the system for x with -y's forcing vector
[y_sim_x_neg_y, t_sim_x_neg_y] = lsim(sys_x, -yReg(L, r), dt*(L-1), xReg(1, 1:r-1)); 

%Recostuct and simulate the system for y with -x's forcing vector
[y_sim_y_neg_x, t_sim_y_neg_x] = lsim(sys_y, -xReg(L, r), dt*(L-1), yReg(1, 1:r-1));


% ---------------PLOTTING X-------------------

% Plot delay-embedded attractor for x
L = 300:length(tspan)-300;
figure;
subplot(2, 2, 1);
set(gcf, 'Name', 'Delay Embedded Attractor X and Error');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Original delay embedding of x
title('Original Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the reconstructed attractor for x
subplot(2, 2, 2);
plot3(y_sim_x(L,1), y_sim_x(L,2), y_sim_x(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);

% Plot the reconstructed attractor for x with 2nd data
subplot(2, 2, 3);
plot3(y_sim_x2(L,1), y_sim_x2(L,2), y_sim_x2(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Delay Embedded Attractor X with 2nd Data');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);

% Compute error over time
error_x = sqrt(sum((V_x(L,1:r-1) - y_sim_x(L,1:r-1)).^2, 2));
error_x2 = sqrt(sum((V_x2(L,1:r-1) - V_x(L,1:r-1)).^2, 2));

% Plot error over time
subplot(2, 2, 4);
plot(tspan(L), error_x, 'r-', 'LineWidth', 1.5);
hold on;
plot(tspan(L), error_x2, 'b-', 'LineWidth', 1.5);
%plot3(y_sim_x2(L2,1),y_sim_x2(L2,2),y_sim_x2(L2,3), 'r-', 'LineWidth', 1.5);
title('Error Between Original and Reconstructed X Over Time');
xlabel('Time');
ylabel('Error');
grid on;
axis tight;


% % ---------------PLOTTING Y-------------------
% % Plot delay-embedded attractor for y
% figure;
% subplot(1, 2, 1);
% set(gcf, 'Name', 'Delay Embedded Attractor');
% set(gcf, 'NumberTitle', 'off');
% plot3(V_y(:,1), V_y(:,2), V_y(:,3));  % Using V_y for the delay embedding of y
% title('Delay Embedded Attractor Y');
% view(-15, 65);

% % Plot the reconstructed attractor for y
% subplot(1, 2, 2);
% set(gcf, 'Name', 'Reconstructed Delay Embedded Attractor Y');
% set(gcf, 'NumberTitle', 'off');

% plot3(y_sim_y(L,1), y_sim_y(L,2), y_sim_y(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
% title('Reconstructed Delay Embedded Attractor Y');
% axis tight;
% xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
% view(-15, 65);

% Compute error over time
error_y = sqrt(sum((V_y(L,1:r-1) - y_sim_y(L,1:r-1)).^2, 2));
error_y2 = sqrt(sum((V_y2(L,1:r-1) - V_y(L,1:r-1)).^2, 2));


% %Plot forcing vector of y before and after the filter
% figure
% % Set the figure's name
% set(gcf, 'Name', 'Forcing Vector of y Before and After the Filter');
% set(gcf, 'NumberTitle', 'off');
% plot(tspan(L), yReg(L,4), 'r', 'LineWidth', 0.5)
% hold on;
% plot(tspan(L), applyKalmanFilter(yReg(L,4)), 'b', 'LineWidth', 0.5)
% legend("Original", "Filtered")
% box on;



% ---------------PLOTTING X WITH Y FORCING VECTOR-------------------
% Plot delay-embedded attractor for x
figure;
subplot(1, 2, 1);
set(gcf, 'Name', 'Delay Embedded Attractor X with Y Forcing and Error');
set(gcf, 'NumberTitle', 'off');
plot3(V_x(:,1), V_x(:,2), V_x(:,3));  % Original delay embedding of x
title('Original Delay Embedded Attractor X');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
view(-15, 65);

% Plot the reconstructed attractor for x with y's forcing vector
subplot(1, 2, 2);
plot3(y_sim_x_y(L,1), y_sim_x_y(L,2), y_sim_x_y(L,3), 'Color', [0 0 0.5], 'LineWidth', 1.5);
title('Reconstructed Attractor X with Y Forcing');
xlabel('v_1'), ylabel('v_2'), zlabel('v_3');
axis tight;
view(-15, 65);


error_x_y = sqrt(sum((V_x(L,1:r-1) - y_sim_x_y(L,1:r-1)).^2, 2));
error_x_neg_y = sqrt(sum((V_x(L,1:r-1) - y_sim_x_neg_y(L,1:r-1)).^2, 2));
%error_x_y = sqrt(sum((V_x(L,1:3) - y_sim_x_y(L,1:3)).^2, 2));

% % Plot error over time
% subplot(1, 3, 3);
% plot(tspan(L), error_x_y, 'r-', 'LineWidth', 1.5);
% title('Error Between Original X and Reconstructed X with Y Forcing');
% xlabel('Time');
% ylabel('Error');
% grid on;
% axis tight;


% -------------COMPARING X vs XY vs X_neg_y ERROR----------------
% plot error x vs error xy on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'Error X vs Error XY vs Error X_neg_y');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L), error_x, 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L), error_x_y, 'b', 'LineWidth', 0.5)
plot(tspan(L), error_x_neg_y, 'm', 'LineWidth', 0.5)
legend("Error X", "Error XY", "Error X(-Y)")
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
error_y_neg_x = sqrt(sum((V_y(L_plot,1:r-1) - y_sim_y_neg_x(L_plot,1:r-1)).^2, 2));

% Plot error over time
subplot(1, 3, 3);
plot(tspan(L_plot), error_y_x, 'r-', 'LineWidth', 1.5);
title('Error Between Original Y and Reconstructed Y with X Forcing');
xlabel('Time');
ylabel('Error');
grid on;
axis tight;

% -------------COMPARING Y vs YX vs Y_neg_x ERROR----------------
% plot error y vs error yx on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'Error Y vs Error YX vs Error Y_neg_x');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), error_y, 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), error_y_x, 'b', 'LineWidth', 0.5)
plot(tspan(L_plot), error_y_neg_x, 'm', 'LineWidth', 0.5)
legend("Error Y", "Error YX", "Error Y(-X)")
box on;

% -------------COMPARING X vs XY and Y vs YX----------------
% plot x vs x2 and y vs y2 on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'X vs XY vs X(-Y)');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_x(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_x_y(L_plot,1), 'b', 'LineWidth', 0.5)
%plot(tspan(L_plot), y_sim_x_neg_y(L_plot,1), 'm', 'LineWidth', 0.5)
legend("X", "XY", "X(-Y)")
box on;

figure
% Set the figure's name
set(gcf, 'Name', 'Y vs YX vs Y(-X)');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_y(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_y_x(L_plot,1), 'b', 'LineWidth', 0.5)
%plot(tspan(L_plot), y_sim_y_neg_x(L_plot,1), 'm', 'LineWidth', 0.5)
legend("Y", "YX", "Y(-X)")
box on;

%----------------COMPARING X AND X2----------------
% plot x vs x2 on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'X vs X2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_x(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_x2(L_plot,1), 'b', 'LineWidth', 0.5)
legend("X", "X2")
box on;

%----------------COMPARING Y AND Y2----------------
% plot y vs y2 on one grpah
figure
% Set the figure's name
set(gcf, 'Name', 'Y vs Y2');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L_plot), y_sim_y(L_plot,1), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L_plot), y_sim_y2(L_plot,1), 'b', 'LineWidth', 0.5)
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


L = 300:length(tspan)-300;
%Plot v11 of y vs v11 of x on 1 graph
figure
% Set the figure's name
set(gcf, 'Name', 'V11y vs V11x');
set(gcf, 'NumberTitle', 'off');
plot(tspan(L), yReg(L,11), 'r', 'LineWidth', 0.5)
hold on
plot(tspan(L), -xReg(L,11), 'b', 'LineWidth', 0.5)
legend("V11y", "-V11x")
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