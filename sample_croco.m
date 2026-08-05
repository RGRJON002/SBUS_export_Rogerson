%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE_CROCO
%
% Sample a 3-D Eulerian CROCO field along Lagrangian particle trajectories.
%
% INPUTS
%   croco_var   : [Ny x Nx x Nt] Eulerian variable
%   croco_lon   : [Ny x Nx]
%   croco_lat   : [Ny x Nx]
%   time_croco  : [Nt x 1] datetime
%
%   float_lon   : [Nfloats x 31 x 31]
%   float_lat   : [Nfloats x 31 x 31]
%   time_float  : [31 x 31] datetime
%
% OUTPUT
%   OUT         : [Nfloats x 30 x 31]
%
% Written: Jonathan Rogerson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OUT = sample_croco( ...
    croco_var, croco_lon, croco_lat, time_croco, ...
    float_lon, float_lat, time_float)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nFloat,~,nRelease] = size(float_lon);

OUT = NaN(nFloat,30,nRelease);

% Flatten CROCO grid once
croco_lon_flat = croco_lon(:);
croco_lat_flat = croco_lat(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-compute temporal indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Pre-computing temporal indices...')

time_croco_day = dateshift(time_croco,'start','day');

time_index = zeros(30,nRelease);

for r = 1:nRelease

    for d = 2:31

        target_day = dateshift(time_float(d,r),'start','day');

        [~,time_index(d-1,r)] = min(abs(time_croco_day - target_day));

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample CROCO field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Sampling Eulerian field...')

for r = 1:nRelease

    fprintf('Release %02d of %02d\n',r,nRelease)

    for d = 2:31

        % Corresponding CROCO time slice
        tind = time_index(d-1,r);

        % Flatten the CROCO field once for this day
        croco_slice = croco_var(:,:,tind);
        croco_flat  = croco_slice(:);

        for f = 1:nFloat

            lon = float_lon(f,d,r);
            lat = float_lat(f,d,r);

            % Ignore missing particles
            if isnan(lon) || isnan(lat)
                continue
            end

            % Distance to every CROCO grid cell
            dist = hypot(croco_lon_flat-lon, ...
                         croco_lat_flat-lat);

            % Nearest neighbour
            [~,ind] = min(dist);

            % Sample variable
            OUT(f,d-1,r) = croco_flat(ind);

        end

    end

end

disp('Sampling complete.')

end