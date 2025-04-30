function [t, sol] = get_lorenz_data(initial_conditions, tspan)
        % Parameters
        sigma = 10;
        rho = 28;
        beta = 8/3;
    
        % Define the Lorenz system as a function
        lorenz = @(t, xyz) [sigma*(xyz(2) - xyz(1)); ...
                            xyz(1)*(rho - xyz(3)) - xyz(2); ...
                            xyz(1)*xyz(2) - beta*xyz(3)];
    
        % Solve using ode45
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        [t, sol] = ode45(lorenz, tspan, initial_conditions, options);
    
    end
    