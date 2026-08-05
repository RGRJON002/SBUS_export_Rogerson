%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute statistics of a Lagrangian sampled variable as a function of age
%
% INPUT
%   VAR : [Nfloat x Nage x Nrelease]
%
% OUTPUT
%   stats.mean
%   stats.std
%   stats.median
%   stats.n
%
% Jonathan Rogerson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats = age_statistics(VAR)

nAge = size(VAR,2);

stats.mean   = NaN(1,nAge);
stats.std    = NaN(1,nAge);
stats.median = NaN(1,nAge);
stats.n      = NaN(1,nAge);

for d = 1:nAge

    % Collect all floats and releases for this age
    tmp = VAR(:,d,:);
    tmp = tmp(:);

    % Remove NaNs (particles older than coastal age)
    tmp = tmp(~isnan(tmp));

    if isempty(tmp)
        continue
    end

    stats.mean(d)   = mean(tmp);
    stats.std(d)    = std(tmp);
    stats.median(d) = median(tmp);
    stats.n(d)      = numel(tmp);

end

end