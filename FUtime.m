%%%% Function to convert 6-hourly trajectory times to daily times.
%%%% Retains the initial release time (Day 0) and averages every
%%%% subsequent group of four 6-hourly outputs.
%%%%
%%%% Input:
%%%%     input : [121 x Fcount] datetime or numeric array
%%%%
%%%% Output:
%%%%     output : [31 x Fcount]

function output = FUtime(input)

    % Get dimensions
    [dim1, dim2] = size(input);

    % Check that the input is compatible
    if mod(dim1 - 1,4) ~= 0
        error('The first dimension minus the first row must be divisible by 4.');
    end

    % Retain the initial release time (Age = 0)
    first_row = input(1,:);

    % Average each subsequent day (4 x 6-hourly outputs)
    reshaped = reshape(input(2:end,:),4,[],dim2);

    % Mean over the 6-hourly samples
    averaged = squeeze(mean(reshaped,1));

    % Concatenate Day 0 with daily means
    output = [first_row; averaged];

end