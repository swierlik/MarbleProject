function y_reconstructed = reconstructHAVOK(A, B, xReg, dt)
    % Ensure consistency in dimensions between B and D
    D = zeros(size(A, 1), size(B, 2)); % Correctly size D
    
    % Create state-space system
    sys = ss(A, B, eye(size(A, 1)), D);
    
    % Time vector
    t_vector = dt * (0:size(xReg, 1)-1);
    
    % Simulate system
    [y_reconstructed, ~] = lsim(sys, xReg(:, end), t_vector, xReg(1, 1:end-1));
end
