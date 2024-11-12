% Lorenz Attractor Animation with Multiple Trajectories
% Parameters
sigma = 10;   % Prandtl number
rho = 28;     % Rayleigh number
beta = 8/3;   % Geometric factor

% Define the Lorenz system as a function
lorenz = @(t, X) [ sigma * (X(2) - X(1));
                   X(1) * (rho - X(3)) - X(2);
                   X(1) * X(2) - beta * X(3)];

% Time span
tspan = [0 50];            % Simulation time

% Number of trajectories
numTrajectories = 10;
initialConditions = zeros(numTrajectories, 3);

% Generate slightly different initial conditions
rng(1); % Seed for reproducibility
for i = 1:numTrajectories
    initialConditions(i, :) = [-8 + 0.1*randn, 8 + 0.1*randn, 27 + 0.1*randn];
end

% Prepare the figure for animation
figure;
hold on;
colors = lines(numTrajectories); % Different colors for each trajectory

% Plotting handles for paths and points
pathHandles = gobjects(numTrajectories, 1);
pointHandles = gobjects(numTrajectories, 1);

for i = 1:numTrajectories
    % Solve the system of equations for each set of initial conditions
    [t, X] = ode45(lorenz, tspan, initialConditions(i, :));
    
    % Initialize path and point handles
    pathHandles(i) = plot3(X(1,1), X(1,2), X(1,3), 'Color', colors(i,:), 'LineWidth', 1.5);
    pointHandles(i) = plot3(X(1,1), X(1,2), X(1,3), 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
end

title('Lorenz Attractor Animation with Multiple Trajectories');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis tight;
view(45, 45); % Adjust the viewing angle for better visualization

% Animation loop
for k = 1:length(t)
    for i = 1:numTrajectories
        % Update the path for each trajectory
        [t, X] = ode45(lorenz, [0 tspan(2)], initialConditions(i, :));
        set(pathHandles(i), 'XData', X(1:k,1), 'YData', X(1:k,2), 'ZData', X(1:k,3));
        
        % Update the moving point for each trajectory
        set(pointHandles(i), 'XData', X(k,1), 'YData', X(k,2), 'ZData', X(k,3));
    end
    
    % Pause to create animation effect
    pause(0.01); % Adjust pause duration for animation speed
end
