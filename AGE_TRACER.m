%%%% This function takes the float lon and lot positions and uses the
%%%% bathymetry depth as an indicator of coastal waters in origin after
%%%% they have been backwards integrated in time. The function takes the
%%%% float_lon and float_lat outputs from ROFF as well as the topography
%%%% data from one of the croco files and the grid dimensions.
%%%% Written: Jonathan Rogerson

%%

function age_matrix = AGE_TRACER(float_lon,float_lat, croco_lon, croco_lat, croco_top)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACER AGE
% To create the age tracer, need to loop over the lon and lat array. 

    % Preallocate depth matrix with NaNs
    depth_matrix = NaN(size(float_lon));

    % Flatten croco_lon and croco_lat for easier searching
    croco_lon_flat = croco_lon(:);
    croco_lat_flat = croco_lat(:);
    croco_top_flat = croco_top(:);

    % Loop through each float coordinate
    for i = 1:size(float_lon, 1)
        for j = 1:size(float_lon, 2)
            for k = 1:size(float_lon, 3)
                % Get float coordinates
                lon = float_lon(i, j, k);
                lat = float_lat(i, j, k);

                % Skip if NaN
                if isnan(lon) || isnan(lat)
                    continue;
                end

                % Compute Euclidean distance to all croco points
                distances = sqrt((croco_lon_flat - lon).^2 + (croco_lat_flat - lat).^2);

                % Find the index of the minimum distance
                [~, min_idx] = min(distances);

                % Assign the corresponding depth value
                depth_matrix(i, j, k) = croco_top_flat(min_idx);
            end
        end
    end

    % Compute the age tracer

    age_matrix = NaN(size(float_lon,1),size(float_lon,3));

    for i = 1:size(float_lon, 1)
        for j = 1:size(float_lon, 3)
            tmp = depth_matrix(i,:,j);  % Single float trajecotry
            tmp(:,1) = [];  % Remove initial index
          % Skip if NaN
            if isnan(tmp(1)) 
                continue;
            end

            id=find(tmp<200,1,'first');

             % Skip if empty 
            if isempty(id) 
                continue;
            end

            age_matrix(i, j) = id;
        end
    end
end