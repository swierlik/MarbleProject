function x_reconstructed = dehankelize(H)
    % Get dimensions of Hankel matrix
    [m, n] = size(H);
    
    % Allocate space for the reconstructed signal
    x_reconstructed = zeros(m + n - 1, 1);
    count = zeros(m + n - 1, 1);

    % Compute the sum along anti-diagonals
    for i = 1:m
        for j = 1:n
            index = i + j - 1;
            x_reconstructed(index) = x_reconstructed(index) + H(i, j);
            count(index) = count(index) + 1;
        end
    end

    % Take the average to correct for overlap
    x_reconstructed = x_reconstructed ./ count;
end
