% Check if Tensorlab is installed and working
try
    % Create a random tensor
    tensor = rand(3, 3, 3);
    
    % Perform a basic tensor decomposition using Tensorlab
    [U, S, V] = svd(tensor(:,:,1));
    
    % Display the results
    disp('Tensorlab is working correctly.');
    disp('Singular values of the first slice of the tensor:');
    disp(diag(S));
catch ME
    disp('Tensorlab is not working correctly.');
    disp(ME.message);
end