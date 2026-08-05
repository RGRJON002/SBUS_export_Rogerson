function output = FU3D(input)
    % Get the size of the input matrix
    [dim1, dim2, dim3] = size(input);
    
    % Ensure the second dimension can be properly averaged
    if mod(dim2 - 1, 4) ~= 0
        error('The second dimension minus the first column must be divisible by 4.');
    end
    
    % Preserve the first column
    first_col = input(:, 1, :);
    
    % Reshape the second dimension (excluding first column) to group every 4 elements
    reshaped = reshape(input(:, 2:end, :), dim1, 4, [], dim3);
    
    % Compute the mean across the second dimension (group of 4)
    averaged = mean(reshaped, 2);
    
    % Squeeze to remove the singleton dimension from averaging
    averaged = squeeze(averaged);
    
    % Concatenate the preserved first column with the averaged data
    output = cat(2, first_col, averaged);
end