%%%% Script to extract depth of a chosen variable/tracer
%%%% Only interested in January. The correspinding period of the Expriments
%%%% From the avg files (not diabio), variables to extract are:
%%%% POC concentration - DET - DET
%%%% Chlorophyll - surface (make depth -1)
%%%% Nurtient load - NO3 - NO3 (from avg file)
%%%% Phytoplankton biomass - PHY - PHY
%%%% units are mmol N m^-3
%%%% This script is a lot simpler than the flux script as we are getting
%%%% a concentration.
%%%% Written: Jonathan Rogerson
%%%% Date: March 2025
%%%% EDIT
%%%% We are going to extract the remineralisation rate at 30 m by using an interpolation. This is correct because we
%%%% are interested in the 
%%%% Remineralisation rate - NFlux_ReminD - REM

addpath('/usr/local/MATLAB/R2025b/toolbox/matlab/imagesci_utils/')

CROCO_path = '/media/jrogerson/JONO/CROCO_BIO';
Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

% SEASON
MS = 7;      % Summer = January, Winter = July

if MS == 1
    months = [12 1];
elseif MS == 7
    months = [6 7];
else
    error('MS must be 1 (summer) or 7 (winter).')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc');
% Read in the mask
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
% Read in topo
h = ncread(file,'h');
h = h.*mask;
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(file,'lat_rho');
CROCO_lon=ncread(file,'lon_rho');

%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%

var = 'NFlux_ReminD';
type = 'r';           % rho grid
depth = -30;          % Depth [m]
var_type = 'REM';     % Output name

start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VAR_mean = [];
TIME = [];

for i = Ymin:Ymax

    %--------------------------------------------------------------
    % Determine which months to read
    % Summer : December (previous year) + January (current year)
    % Winter : June + July (current year)
    %--------------------------------------------------------------

    if MS == 1
        years_to_read  = [i-1, i];
        months_to_read = [12, 1];
    elseif MS == 7
        years_to_read  = [i, i];
        months_to_read = [6, 7];
    else
        error('MS must be 1 (summer) or 7 (winter).')
    end

    %--------------------------------------------------------------
    % Read the required months
    %--------------------------------------------------------------

    w_slice = [];

    for m = 1:2

        file = strcat(CROCO_path,'/','diabio_avg_Y', ...
            string(years_to_read(m)),'M',string(months_to_read(m)),'.nc');

        if exist(file,'file')

            disp(file)
            disp('Reading data')
            tmp_slice = [];
            
            %--------------------------------------------------------------
            % Read time
            % If this is a daibio file, read time from the corresponding avg file 
            %--------------------------------------------------------------
            time_file = file;

            if contains(time_file,'diabio')
                time_file = strrep(time_file,'diabio_','');
            end
            
            time = ncread(time_file,'time');
            
            for k = 1:length(time)
            
                tmp_slice(:,:,k) = get_hslice( ...
                    convertStringsToChars(file), ...
                    convertStringsToChars(file), ...
                    var,k,depth,type);
            
            end

            % Store
            disp('Storing variable and time')

            TIME = cat(1,TIME,time);
            w_slice = cat(3,w_slice,tmp_slice);

            fprintf('%s\n', file)
            fprintf('length(time) = %d\n', length(time))
            fprintf('size(w_slice,3) = %d\n\n', size(tmp_slice,3))

        else

            warning(['Missing file: ',file])

        end

    end

    % Store this year's Dec-Jan (or Jun-Jul) cube
    VAR_mean = cat(3,VAR_mean,w_slice);

end

% Format time

% [~, ia, ~] = unique(TIME);
% VAR_mean = VAR_mean(:,:,ia);
% TIME = TIME(ia);

% CHECK

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save yearly trajectory fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = Ymin:Ymax;

for y = years

    if MS == 1
        % Summer: December (previous year) + January (current year)
        idx = (year(CROCO_time) == y-1 & month(CROCO_time) == 12) | ...
              (year(CROCO_time) == y   & month(CROCO_time) == 1);

    elseif MS == 7
        % Winter: June + July (current year)
        idx = (year(CROCO_time) == y & month(CROCO_time) == 6) | ...
              (year(CROCO_time) == y & month(CROCO_time) == 7);

    end

    % Extract this trajectory year's data
    VAR = VAR_mean(:,:,idx);
    time_year = CROCO_time(idx);

    % Save
    filename = strcat(var_type,'_',string(y),'.mat');
    save(filename,'VAR','time_year','-v7.3');

    disp(['Saved ',filename])

end
