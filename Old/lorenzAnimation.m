% Load the necessary data
load('Data/lorenzData.mat');          % Contains 'sol', 't', 'dt'
load('Data/systemData.mat');          % Contains 'xReg', 'r', 'tspan', 'dt'

% Extract x, y, z trajectories and forcing vector
x = sol(:, 1);                        % Extract x from the system
y = sol(:, 2);                        % Extract y from the system
z = sol(:, 3);                        % Extract z from the system
forcing_vector = xReg(:, r);          % Forcing vector

% Define the time vector
L = 300:(length(xReg) - 300);         % Exclude transients
t = tspan(L);                         % Corresponding time vector

% Prepare for the animation
figure;

% Create subplots for system trajectory and forcing vector
subplot(1, 2, 1);
h1 = plot3(x(1), y(1), z(1), 'b', 'LineWidth', 1.5); % Trajectory line
hold on;
h2 = plot3(x(1), y(1), z(1), 'ro', 'MarkerFaceColor', 'r'); % Current point
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('System Trajectory');
view(3);

subplot(1, 2, 2);
h3 = plot(t(1), forcing_vector(1), 'r', 'LineWidth', 1.5); % Forcing line
hold on;
h4 = plot(t(1), forcing_vector(1), 'ko', 'MarkerFaceColor', 'k'); % Current forcing
xlabel('Time');
ylabel('Forcing Magnitude');
grid on;
title('Forcing Vector');

% Animation loop
for k = 2:length(L)
    % Update trajectory plot
    set(h1, 'XData', x(1:k), 'YData', y(1:k), 'ZData', z(1:k));
    set(h2, 'XData', x(k), 'YData', y(k), 'ZData', z(k));
    
    % Update forcing vector plot
    set(h3, 'XData', t(1:k), 'YData', forcing_vector(1:k));
    set(h4, 'XData', t(k), 'YData', forcing_vector(k));
    
    % Pause for smooth animation
    pause(0.01);
end
