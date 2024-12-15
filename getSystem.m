function [V, A_x, B_x, xReg, t_reg] = getSystem(x, hankelRows, r, dt, tspan)
    % Generate Hankel matrix
    HankelMatrix_x = hankel(x(1:hankelRows), x(hankelRows:end));
    [~, ~, V_full] = svd(HankelMatrix_x, 'econ');

    % Check if V_full is empty or ill-conditioned
    if isempty(V_full)
        error('SVD failed: V_full is empty.');
    end

    % Compute derivatives
    num_V = size(V_full, 1);
    dV = zeros(num_V - 5, r);
    for i = 3:num_V - 3
        for k = 1:r
            dV(i-2, k) = (1 / (12 * dt)) * (-V_full(i+2, k) + 8 * V_full(i+1, k) - 8 * V_full(i-1, k) + V_full(i-2, k));
        end
    end

    % Reduced variables from SVD
    xReg = V_full(3:end-3, 1:r);

    % Adjust time vector to match xReg
    t_reg = tspan(hankelRows + 3:end - 2);

    % Ensure lengths match
    if length(t_reg) ~= size(xReg, 1)
        error('Length of t_reg does not match size of xReg.');
    end

    % Build library of nonlinear time series
    polyorder = 1;
    usesine = 0;
    Theta = poolData(xReg, r, polyorder, usesine);

    % Normalize columns of Theta
    normTheta = zeros(size(Theta, 2), 1);
    for k = 1:size(Theta, 2)
        normTheta(k) = norm(Theta(:, k));
        if normTheta(k) ~= 0
            Theta(:, k) = Theta(:, k) / normTheta(k);
        else
            Theta(:, k) = 0;
        end
    end

    % Sparse regression to obtain Xi
    lambda = 0.0;
    Xi = zeros(size(Theta, 2), r - 1);
    for k = 1:r-1
        Xi(:, k) = sparsifyDynamics(Theta, dV(:, k), lambda * k, 1);
    end

    % De-normalize Xi
    for k = 1:length(Xi)
        if normTheta(k) ~= 0
            Xi(k, :) = Xi(k, :) / normTheta(k);
        else
            Xi(k, :) = 0;
        end
    end

    % Extract system matrices A and B
    if size(Xi, 1) < r+1
        error('Xi does not have enough rows.');
    end
    A_x_full = Xi(2:r+1, 1:r-1)';
    if size(A_x_full, 2) < r-1
        error('A_x_full does not have enough columns.');
    end
    B_x = A_x_full(:, end);
    A_x = A_x_full(:, 1:end-1);

    % Return the first r columns of V_full for consistency
    V = V_full(:, 1:r);
end
