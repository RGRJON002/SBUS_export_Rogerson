%%%% Script to extract depth of a chosen variable
%%%% For this work, we want to extract new production and export of POC as
%%%% indices from the January files. Also want to calculate open ocen
%%%% Compute open ocean values of chlrophyll and POC

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

myCHL = [];
DET = [];

for i = Ymin:Ymax
    % DECLARE the files to loop over
    % CROCO file
    CROCO_file = strcat('avg_Y',string(i),'M1.nc');
    % DET file
    DET_file = strcat('DET_',string(i),'.mat');
    % CROCO CHL
    file = strcat(CROCO_path,'/',CROCO_file);
    % Read in the mask
    mask=ncread(file,'mask_rho');
    mask(mask==0)=nan;
    % Read in topo
    h = ncread(file,'h');
    h = h.*mask;
    % Define my region, will focus on the coastal domain
    CROCO_lat=ncread(file,'lat_rho');
    CROCO_lon=ncread(file,'lon_rho');
    % CHL
    CHL = double(ncread(file,'CHLA'));
    CHL = squeeze(CHL(:,:,50,:));
    CHL = CHL.*mask;
    myCHL = cat(3,myCHL,CHL);
    
    % DET
    load(strcat('DET_',string(i),'.mat'));
    det = VAR;
    det = permute(det, [2, 1, 3]);
    det = det.*mask;
    DET = cat(3,DET,det);
end
    
%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%

% Calculate NaN
myCHL = mean(myCHL,3);
DET = mean(DET,3);

% Convert units
DET = DET * 6.625;

mask = h' > 1000;

% Extract chlorophyll values where depth > 1000
det_filtered = DET(mask);
chl_filtered = myCHL(mask);

% Compute the mean, ignoring NaNs if present
det_filtered = mean(det_filtered, 'omitnan');
chl_filtered = mean(chl_filtered, 'omitnan');


% Display result
disp(['Mean Chlorophyll where depth > 1000 m: ', num2str(det_sink)])



