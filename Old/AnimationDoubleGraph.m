% Lorenz Attractor Animation with Trailing Path and X vs Time Plot
% Parameters
sigma = 10;   % Prandtl number
rho = 28;     % Rayleigh number
beta = 8/3;   % Geometric factor

% Define the Lorenz system as a function
lorenz = @(t, X) [ sigma * (X(2) - X(1));
                   X(1) * (rho - X(3)) - X(2);
                   X(1) * X(2) - beta * X(3)];

% Time span and initial conditions
tspan = [0 50];            % Simulation time
X0 = [1; 1; 1];            % Initial conditions [x0, y0, z0]

% Solve the system of equations using ode45
[t, X] = ode45(lorenz, tspan, X0);

% Set up the figure with subplots
figure;

% 3D Trajectory Plot
subplot(1, 2, 1);
h = plot3(X(1,1), X(1,2), X(1,3), 'b', 'LineWidth', 1.5); % Initial path
hold on;
p = plot3(X(1,1), X(1,2), X(1,3), 'ro', 'MarkerFaceColor', 'r'); % Moving point
title('Lorenz Attractor');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis tight;
view(45, 45); % Adjust the viewing angle for better visualization

% X vs Time Plot
subplot(1, 2, 2);
plot(t, X(:,1), 'b', 'LineWidth', 1.5);
title('X vs Time');
xlabel('Time');
ylabel('X');
grid on;

% Animation loop
for k = 1:length(t)
    % Update the 3D path
    set(h, 'XData', X(1:k,1), 'YData', X(1:k,2), 'ZData', X(1:k,3));
    
    % Update the moving point
    set(p, 'XData', X(k,1), 'YData', X(k,2), 'ZData', X(k,3));
    
    % Update the X vs Time plot
    subplot(1, 2, 2);
    plot(t, X(:,1), 'b', 'LineWidth', 1.5); % Full trajectory
    hold on;
    plot(t(k), X(k,1), 'ro'); % Current point
    title('X vs Time');
    xlabel('Time');
    ylabel('X');
    grid on;
    
    % Pause to create animation effect
    subplot(1, 2, 1); % Switch back to 3D plot
    pause(0.01); % Adjust pause duration for animation speed
end


