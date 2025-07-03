%%%% For a particular variable, it would be nice to co-index it with
%%%% specific age/age range waters. So this function allows a variable (2D)
%%%% to be decalred as well as a declaraetion of waters only older than a
%%%% chosen value to then be coextracted.
%%%% Written: Jonathan Rogerson

function [myVAR, age_vals, lon_vals, lat_vals] = agetovar(float_lon, float_lat, age_array,myage, myvar, croco_lon, croco_lat) 

% Flatten croco_lon and croco_lat for easier searching
croco_lon_flat = croco_lon(:);
croco_lat_flat = croco_lat(:);

% Overlay
% Apply a basic filter to eclude very recent water
age_array(age_array<=myage) = NaN;

% Remove NaN values to only plot what I need
% Remove NaN values (if any exist)
lon_vals = squeeze(float_lon(:,1,1));
lat_vals = squeeze(float_lat(:,1,1));
age_vals = age_array;

valid_idx = ~isnan(age_vals);
lon_vals = lon_vals(valid_idx);
lat_vals = lat_vals(valid_idx);
age_vals = age_vals(valid_idx);

% Small loop

myvar_flat = myvar(:);

myVAR = NaN(size(lon_vals));
for i = 1:length(lon_vals)
        % Compute Euclidean distance to all croco points
        distances = sqrt((croco_lon_flat - lon_vals(i)).^2 + (croco_lat_flat - lat_vals(i)).^2);
        % Find the index of the minimum distance
        [~, min_idx] = min(distances);
        % Assign the corresponding depth value
        myVAR(i) = myvar_flat(min_idx);
end

end
