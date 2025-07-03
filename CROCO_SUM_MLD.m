%%%% Compute the mean MLD depth for January over 2011-2015

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the CROCO data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING  GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define my region, will focus on the coastal domain
croco_top = ncread(strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc'),'h');
croco_top = croco_top.*mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the closest points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon_tmp = CROCO_lon(:,1);
lat_tmp = CROCO_lat(1,:);

Latitude = [-32.30,-32.30];
Longitude = [17.10,18.31]; 

% Starting point

[ind_start_i,ind_start_j,km] = mll2grid(Latitude(end),Longitude(end),CROCO_lat,CROCO_lon);

% End point

[ind_end_i,ind_end_j,km] = mll2grid(Latitude(1),Longitude(1),CROCO_lat,CROCO_lon);

lonsec = [CROCO_lon(ind_start_j,ind_start_i), CROCO_lon(ind_end_j,ind_end_i)];
lonsec = double(lonsec(1,[1,end]));

latsec = [CROCO_lat(ind_start_j,ind_start_i), CROCO_lat(ind_end_j,ind_end_i)];
latsec = double(latsec(1,[1,end]));
%% a vertical section

start

SALT = [];
TEMP = [];
X = [];
Z =[];
TIME = [];

for i = Ymin:Ymax
    for j = 1
    file = strcat(CROCO_path,'/','avg_Y',string(i),'M',string(j),'.nc');
    file = convertStringsToChars(file);
    if exist(file, 'file')
        disp(file)
        disp('Reading data')
        % Get the time domain 
        time = ncread(file,'time');
        % Get temperature
        for k = 1:length(time)
            [x(:,:,k),z(:,:,k),temp(:,:,k)] = get_section(file,file,lonsec,latsec,'temp',k);
            [~,~,salt(:,:,k)] = get_section(file,file,lonsec,latsec,'salt',k);
        end
        %  Store arrays
        TEMP = cat(3,TEMP,temp);
        SALT = cat(3,SALT,salt);
        X = cat(3,X,x);
        Z = cat(3,Z,z);
        TIME = cat(1,TIME,time);
    else
        disp(strcat('No data for',file))
    end
    clear salt temp time
    end
end

% Quick clean

[~, ia, ~] = unique(TIME);
TEMP = TEMP(:,:,ia);
X = X(:,:,ia);
Z = Z(:,:,ia);
SALT = SALT(:,:,ia);
TIME = TIME(ia);

%% LOAD in DATA
%MODEL_SHBML = struct('X',X,'Z',Z,'SALT',SALT,'TEMP',TEMP,...
%    'TIME',TIME);

%save('MODEL_SHBML.mat','MODEL_SHBML');

load 'MODEL_SHBML.mat'

X = MODEL_SHBML.X;
Z = MODEL_SHBML.Z;
TEMP = MODEL_SHBML.TEMP;
SALT = MODEL_SHBML.SALT;
TIME = MODEL_SHBML.TIME;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE DENSITY and MLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(TIME)
    for j = 1:size(TEMP,2)
        z=flipud(squeeze(Z(:,j,i)));
        t=flipud(squeeze(TEMP(:,j,i)));
        [MLD_T(j,i),qe,imf]=get_mld(z,t);
    end
end

% Compute mean profile

MLD_mean = mean(MLD_T, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME index the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);