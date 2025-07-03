%%%% Script to extract depth of a chosen variable/tracer
%%%% Only interested in January. The correspinding period of the Expriments
%%%% From the avg files (not diabio), variables to extract are:
%%%% POC concentration - DET - DET
%%%% Zooplankton biomass - ZOO - ZOO
%%%% units are mmol N m^-3
%%%% This script is a lot simpler than the flux script as we are getting
%%%% a concentration.
%%%% Written: Jonathan Rogerson
%%%% Date: March 2025


addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

CROCO_path = '/media/jono/JONO/CROCO_BIO';
Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

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

var = 'ZOO';
type = 'r'; % rho grid
depth = -30;  % Depth in [m]
var_type = 'ZOO'; % output name

start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VAR_mean = [];
TIME = [];
for i = Ymin:Ymax
    file = strcat(CROCO_path,'/','avg_Y',string(i),'M',string(1),'.nc');
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Read in the time array
        time = ncread(strcat(CROCO_path,'/','avg_Y',string(i),'M',string(1),'.nc'),'time');
        % Read in the vertical velocity data
        addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')
        w = ncread(file,var);
        for k = 1:length(time)
            w_slice(:,:,k) = get_hslice(convertStringsToChars(file),...
                convertStringsToChars(file),var,k,depth,type);
        end
        %  Store arrays
        disp('Storing w and time')
        TIME = cat(1,TIME,time);
        VAR_mean = cat(3,VAR_mean,w_slice);
    else
        disp(strcat('No data for',file))
    end
    clear w w_slice
end

% Format time

[~, ia, ~] = unique(TIME);
VAR_mean = VAR_mean(:,:,ia);
TIME = TIME(ia);

% CHECK

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);

% Save to a structure

% Save to a structure
years = Ymin:Ymax;

for y = years
    % Logical index for current year
    idx = year(CROCO_time) == y;

    % Subset the data
    VAR = VAR_mean(:,:,idx);
    time_year = CROCO_time(idx);

    % Save to a .mat file named by the year
    filename = strcat(var_type,'_',string(y),'.mat');
    save(filename, 'VAR', 'time_year', '-v7.3');  % -v7.3 is safer for large files
end