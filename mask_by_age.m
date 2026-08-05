function VAR = mask_by_age(VAR, coastal_age)

%----------------------------------------------------------------------
% Converts backward trajectories into coastal-age trajectories.
%
% Input:
%   VAR          : [nFloat x nDay x nRelease]
%   coastal_age  : [nFloat x nRelease]
%
% Output:
%   VAR          : same size
%
% After processing:
%
%   VAR(:,1,:)  = first day at the coast
%   VAR(:,2,:)  = one day after leaving the coast
%   ...
%
%----------------------------------------------------------------------

[nFloat,nDay,nRelease] = size(VAR);

OUT = NaN(size(VAR));

for r = 1:nRelease
    for f = 1:nFloat

        age = coastal_age(f,r);

        if isnan(age)
            continue
        end

        % Skip impossible ages
        if age > nDay
            continue
        end

        % Keep only trajectory up to the coast
        traj = VAR(f,1:age,r);

        % Reverse so age starts at the coast
        traj = fliplr(traj);

        % Store beginning at Day 1
        OUT(f,1:length(traj),r) = traj;

    end
end

VAR = OUT;

end
