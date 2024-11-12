% Lorenz Attractor Animation with Trailing Path and Y vs Time Plot
% Load the previously generated Lorenz data
load('lorenzData.mat');

% Downsample factor (reduce the number of points)
downsample_factor = 50;  % Choose an appropriate factor (e.g., 10, 20, etc.)

% Downsample the time and solution data
t_downsampled = t(1:downsample_factor:end);
x_downsampled = sol(1:downsample_factor:end, 1); % Extracting x
y_downsampled = sol(1:downsample_factor:end, 2); % Extracting y
z_downsampled = sol(1:downsample_factor:end, 3); % Extracting z

% Set up the figure with subplots
figure;

% 3D Lorenz Attractor Plot (X, Y, Z)
subplot(1, 2, 1);
h = plot3(x_downsampled(1), y_downsampled(1), z_downsampled(1), 'b', 'LineWidth', 1.5); % Initial path
hold on;
p = plot3(x_downsampled(1), y_downsampled(1), z_downsampled(1), 'ro', 'MarkerFaceColor', 'r'); % Moving point
title('Lorenz Attractor');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis tight;
view(45, 45); % Adjust the viewing angle for better visualization

% Y vs Time Plot
subplot(1, 2, 2);
plot(t_downsampled, y_downsampled, 'b', 'LineWidth', 1.5);
title('Y vs Time');
xlabel('Time');
ylabel('Y');
grid on;

% Animation loop for visualizing the path and Y vs Time evolution
for k = 1:length(t_downsampled)
    % Update the 3D path (x, y, z)
    set(h, 'XData', x_downsampled(1:k), 'YData', y_downsampled(1:k), 'ZData', z_downsampled(1:k));
    
    % Update the moving point
    set(p, 'XData', x_downsampled(k), 'YData', y_downsampled(k), 'ZData', z_downsampled(k));
    
    % Update the Y vs Time plot
    subplot(1, 2, 2);
    plot(t_downsampled, y_downsampled, 'b', 'LineWidth', 1.5); % Full trajectory
    hold on;
    plot(t_downsampled(k), y_downsampled(k), 'ro'); % Current point
    title('Y vs Time');
    xlabel('Time');
    ylabel('Y');
    grid on;
    
    % Pause to create animation effect
    subplot(1, 2, 1); % Switch back to 3D plot
    pause(0.01); % Adjust pause duration for animation speed
end
