%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGE_TRACER
%
% Computes the coastal age of each particle using the sampled bathymetry.
%
% INPUT
%   DPT : [Nfloats x 30 x Nrelease]
%         Lagrangian sampled bathymetry (m)
%
% OUTPUT
%   coastal_age : [Nfloats x Nrelease]
%
%                 Number of days required for each particle to first
%                 encounter water shallower than 200 m.
%
%                 NaN indicates the particle never reaches the shelf.
%
% Written: Jonathan Rogerson
% NOTE: This builds off an earlier version of the code (used to make Fig.
% 5). We just adapt the code here to leverage the bathymtry that is
% actually sampled along the lagrngian trajecotries rather than having to
% recompute a depth matrix for sampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coastal_age = AGE_TRACER(DPT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nFloat,nDay,nRelease] = size(DPT);

coastal_age = NaN(nFloat,nRelease);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute coastal age
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:nRelease

    for f = 1:nFloat

        % Bathymetry history for one particle
        depth = squeeze(DPT(f,:,r));

        % Skip incomplete trajectories
        if all(isnan(depth))
            continue
        end

        % First day entering the shelf (<200 m)
        idx = find(depth < 200,1,'first');

        if ~isempty(idx)
            coastal_age(f,r) = idx;
        end

    end

end

end