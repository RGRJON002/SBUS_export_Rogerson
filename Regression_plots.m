%%%% The objective of this script is to pull together the age computation
%%%% as well as the metrics for NSINK and NPROD and construct three
%%%% regression plots.

%%%% Written: Jonathan Rogerson

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/Documents/PhD/CHAPTER_2/DATA
addpath /home/jono/Documents/PhD/CHAPTER_2/SCRIPTS
addpath /home/jono/Documents/MATLAB/ezyfit2.44/ezyfit
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')

%% MAIN LOOP

% 1) File looping
% 2) Computation
% 3) Data storing

float_path = '/home/jono/Documents/POSTDOC/DATA/Roff_OUT';
CROCO_path = '/media/jono/JONO/CROCO_BIO'; 

Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

Fcount = 31;

% NB

% My age tracer to exclude

cut_age = 0;

% Want to save and store in a stacked array the following processed
% variables
NSINK = [];  % 1-D
NPROD = [];  % 1 -D
myCHL = [];    % 1-D
DET = [];    % 1-D
NREM = [];   % 1-D
NZOO = [];   % 1-D
NGRAZ = [];  % 1-D
AGE = [];    % 1-D
NEXP = [];   % 1-D
AGE_FULL=[]; % 3-D

for i = Ymin:Ymax
    % DECLARE the files to loop over
    %float file
    float_file = strcat('SBUS_floats_bio_Jan',string(i),'.nc');
    % CROCO file
    CROCO_file = strcat('avg_Y',string(i),'M1.nc');
    % NPROD file
    PROD_file = strcat('NPROD_',string(i),'.mat');
    % NSINK file
    SINK = strcat('NSINK_',string(i),'.mat');
    % DET file
    DET_file = strcat('DET_',string(i),'.mat');
    % REM file 
    REM_file = strcat('REM_',string(i),'.mat');
    % NZOO file
    ZOO_file = strcat('NZOO_',string(i),'.mat');
    % NGRAZ file
    GRAZ_file = strcat('NGRAZ_',string(i),'.mat');
    % NEXP file
    NEXP_file = strcat('NEXP_',string(i),'.mat');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Read in the grid and time dimension from the float file

    float_lon = Fto3D(ncread(strcat(float_path,'/',float_file),'lon'),Fcount);
    float_lat = Fto3D(ncread(strcat(float_path,'/',float_file),'lat'),Fcount);
    float_depth = Fto3D(ncread(strcat(float_path,'/',float_file),'depth'),Fcount);
    float_temp = Fto3D(ncread(strcat(float_path,'/',float_file),'temp'),Fcount);
    float_time = ncread(strcat(float_path,'/',float_file),'scrum_time');

    % Process time if it needs to be stacked
    float_time = Ftime(float_time,6674,Fcount);

    % Convert the time array
    time_float = datetime(2000,1,1) + seconds(float_time);

    % Days since initialisation
    % We can be fancy and compute it using date formats but knowing the
    % freq of the saved variable and length will give and indication.

    Duration = (size(float_lat,2) -1)/(24/6); % Number of days floats are tracked 
    freq = 24/6;
    Days = repmat(1:Duration,[freq 1]);
    Days = Days(:)';

    % Pad 
    Days = [Days(1),Days];

    % Get DAYS as same size as other variables

    Days = repmat(Days,size(float_lat,1),1);
    Days = repmat(Days,1,1,size(float_lat,3));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLEAN DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mis_var = 1.0e15;
    ind_mis = find(float_lat == mis_var);
    % Make NaN
    float_lon(ind_mis) = NaN;
    float_lat(ind_mis)= NaN;
    float_depth(ind_mis) = NaN;
    float_temp(ind_mis) = NaN;
    Days(ind_mis) = NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % T_avg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float_lon = FU3D(float_lon);
    float_lat= FU3D(float_lat);
    float_depth = FU3D(float_depth);
    float_temp = FU3D(float_temp);
    Days = FU3D(Days);

    %%%%%%%%%%%%%%%
    % CROCO FILE
    %%%%%%%%%%%%%%

    % Read in the grid dimension of one of the CROCO_files

    croco_lon = double(ncread(strcat(CROCO_path,'/',CROCO_file),'lon_rho'));
    croco_lat = double(ncread(strcat(CROCO_path,'/',CROCO_file),'lat_rho'));
    
    time = ncread(strcat(CROCO_path,'/',CROCO_file),'time');
    CROCO_time = datetime(2000,1,1) + seconds(time);
    % Read in the topography data

    croco_top = ncread(strcat(CROCO_path,'/',CROCO_file),'h');
    mask = ncread(strcat(CROCO_path,'/',CROCO_file),'mask_rho');
    mask(mask==0)=nan;               % Land is = 0, make nan
    croco_top=croco_top.*mask;
    
    % Process the files: CHL, NSINK, NPROD
    
    % CHL 
    CHL = double(ncread(strcat(CROCO_path,'/',CROCO_file),'CHLA'));
    CHL = squeeze(CHL(:,:,50,:));
    CHL = CHL.*mask;
    
    % NSINK
    load(strcat('NSINK_',string(i),'.mat'));
    nsink = VAR;
    nsink = permute(nsink, [2, 1, 3]);
    nsink = nsink.*mask;
    
    % NPROD
    load(strcat('NPROD_',string(i),'.mat'));
    nprod = VAR;
    nprod = permute(nprod, [2, 1, 3]);
    nprod = nprod.*mask;
    
    % DET
    load(strcat('DET_',string(i),'.mat'));
    det = VAR;
    det = permute(det, [2, 1, 3]);
    det = det.*mask;
    
    % NREM
    load(strcat('REM_',string(i),'.mat'));
    rem = VAR;
    rem = permute(rem, [2, 1, 3]);
    rem = rem.*mask;
    
    % NZOO
    load(strcat('NZOO_',string(i),'.mat'));
    zoo = VAR;
    zoo = permute(zoo, [2, 1, 3]);
    zoo = zoo.*mask;
    
    % NGRAZ
    load(strcat('NGRAZ_',string(i),'.mat'));
    graz = VAR;
    graz = permute(graz, [2, 1, 3]);
    graz = graz.*mask;
    
    % NEXP
    load(strcat('NEXP_',string(i),'.mat'));
    exp = VAR;
    exp = permute(exp, [2, 1, 3]);
    exp = exp.*mask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the age
    
    age_matrix = AGE_TRACER(float_lon,float_lat, croco_lon, croco_lat, croco_top);
    
    AGE_FULL = cat(3,AGE_FULL,age_matrix);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find corresponding variables for particular ages of particles
    
    % CHL
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, age_vals, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, CHL(:,:,ind_chl), croco_lon, croco_lat);
        myCHL = cat(1,myCHL,myVAR);
        AGE = cat(1,AGE,age_vals);   % Only needs to be done once as will be identical for the others
    end

    clear myVAR
    % PROD
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, nprod(:,:,ind_chl), croco_lon, croco_lat);
        NPROD = cat(1,NPROD,myVAR);
    end

    clear myVAR
    % SINK
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, nsink(:,:,ind_chl), croco_lon, croco_lat);
        NSINK = cat(1,NSINK,myVAR);
    end
    
    clear myVAR
    % DET
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, det(:,:,ind_chl), croco_lon, croco_lat);
        DET = cat(1,DET,myVAR);
    end
    
    clear myVAR
    % NREM
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, rem(:,:,ind_chl), croco_lon, croco_lat);
        NREM = cat(1,NREM,myVAR);
    end
    
    clear myVAR
    % NZOO
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, zoo(:,:,ind_chl), croco_lon, croco_lat);
        NZOO = cat(1,NZOO,myVAR);
    end
    
    clear myVAR
    % NGRAZ
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, graz(:,:,ind_chl), croco_lon, croco_lat);
        NGRAZ = cat(1,NGRAZ,myVAR);
    end
    
    clear myVAR
    % NEXP
    for j = 1:31
        ind_float = j;
        ind_chl = ind_float + 31 - (2*ind_float) + 1; 
        [myVAR, ~, ~, ~]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),cut_age, exp(:,:,ind_chl), croco_lon, croco_lat);
        NEXP = cat(1,NEXP,myVAR);
    end
    
    clear myVAR
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(strcat('DONE with year --', string(i)));
    
end

%% Regression figure

% Convert units to per day
NPROD = NPROD * 86400;
NSINK = NSINK * 86400;
NREM = NREM * 86400;
NGRAZ = NGRAZ * 86400;
NZOO = NZOO * 86400;
NEXP = NEXP * 86400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try another fit

mean_NPROD = nan(1, 30); % 26 because 31-6+1 = 26 unique values
mean_NSINK = nan(1, 30);
mean_CHL = nan(1, 30);
mean_DET = nan(1,30);
mean_NREM = nan(1,30);
mean_NZOO = nan(1,30);
mean_NGRAZ = nan(1,30);
mean_NEXP = nan(1,30);
PROB = nan(1,30);

for val = 1:30
    % Find indices where array equals the current value
    idx = (AGE == val);
    
    PROB(val) = (sum(double(idx))/length(AGE))*100;
    
    % Calculate means for both arrays if there are valid values
    if any(idx)
        mean_CHL(val)   = mean(myCHL(idx));
        mean_NPROD(val) = mean(NPROD(idx));
        mean_NSINK(val) = mean(NSINK(idx));
        mean_DET(val) = mean(DET(idx));
        mean_NREM(val) = mean(NREM(idx));
        mean_NZOO(val) = mean(NZOO(idx));
        mean_NGRAZ(val) = mean(NGRAZ(idx));
        mean_NEXP(val) = mean(NEXP(idx));
    end
end

%% Basic outline

% All the rates and variables together

xx = 1:30;
names = {'CHL','NSINK','DET','NEXP','NZOO','NREM'};
data = [mean_CHL;mean_NSINK.*-1;mean_DET;mean_NEXP.*-1;mean_NZOO;mean_NREM];

figure
for k = 1:length(names)
    subplot(3,2,k)
    plot(xx,data(k,:),'k--');
    title(names{k})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CtoN = 6.625; % Conversion of carbon:nitorgen according to redfield

% Convert the units to that of mmol C m-2 day-1

mean_NSINK = mean_NSINK*CtoN;
mean_DET = mean_DET*CtoN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nice figure with fitted curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blub = 1:30; % Declare

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1) % CHL
plot(blub,mean_CHL,'k--','LineWidth', 3)
hold on
%pp = 0.79764*exp(-0.13494 * blub) + 1.0214;
%plot(blub, pp, '--', 'Color', 'r', 'LineWidth', 3);
%hold on
%plot(NaN, NaN);
grid on
%lgd = legend('Data',"0.79 exp(-0.13 t) + 1.02",strcat('R^2=',string(0.95)),'Location','best');
%lgd.FontSize = 11.5;
ax = gca;          % Get current axis
ax.XAxis.FontSize = 12;   % Set x-axis tick label size
ax.YAxis.FontSize = 12;   % Set y-axis tick label size
xlabel('Age [Days]','FontSize',13,'FontWeight','bold')
ylabel('Chlorophyll [mg Chla m^{-3}]','FontSize',13,'FontWeight','bold')
title('(a)','FontSize',15,'FontWeight','bold')
xlim([0 30]);

subplot(2,2,2) % DET
plot(blub,mean_DET','k--','LineWidth', 3)
hold on
grid on
ax = gca;          % Get current axis
ax.XAxis.FontSize = 12;   % Set x-axis tick label size
ax.YAxis.FontSize = 12;   % Set y-axis tick label size
xlabel('Age [Days]','FontSize',13,'FontWeight','bold')
ylabel('POC [mmol C m^{-3}]','FontSize',13,'FontWeight','bold')
title('(b)','FontSize',15,'FontWeight','bold')
xlim([0 30]);

subplot(2,2,3)  % NSINK
plot(blub,mean_NSINK*-1,'k--','LineWidth', 3);
hold on
%pp =27.974+ 1.6028*blub - 0.075461*blub.^2;
%plot(blub, pp, '--', 'Color', 'r', 'LineWidth', 3);
%hold on
%plot(NaN, NaN);
grid on
%lgd = legend('Data',"27.97 + 1.60t - 0.08t^2",strcat('R^2=',string(0.84)),'Location','best');
%lgd.FontSize = 11.5;
ax = gca;          % Get current axis
ax.XAxis.FontSize = 12;   % Set x-axis tick label size
ax.YAxis.FontSize = 12;   % Set y-axis tick label size
xlabel('Age [Days]','FontSize',13,'FontWeight','bold')
ylabel('POC vertical flux [mmol C m^{-2} day^{-1}]','FontSize',13,'FontWeight','bold')
title('(c)','FontSize',15,'FontWeight','bold')
xlim([0 30]);

subplot(2,2,4)
plot(blub, mean_NREM','--k','LineWidth',3);
ax = gca;          % Get current axis
ax.XAxis.FontSize = 12;   % Set x-axis tick label size
ax.YAxis.FontSize = 12;   % Set y-axis tick label size
xlabel('Age [Days]','FontSize',13,'FontWeight','bold')
ylabel('Remineralisation rate [mmol C m^{-2} day^{-1}]','FontSize',13,'FontWeight','bold')
title('(d)','FontSize',15,'FontWeight','bold')
grid on
xlim([0 30]);

%%
set(gcf, 'InvertHardcopy', 'off')
print('-f1','Fits_plots','-dpng','-r600');

%% Question 1:

% What is the mean depth (isobath) of floats > 20 days old and north of 35

% Step 1: Identify particles with AGE >= 20 and their depth (isobath)

lon_vals = squeeze(float_lon(:,1,1));
lat_vals = squeeze(float_lat(:,1,1));

lat_vals = repmat(lat_vals,[1, size(AGE_FULL,2), size(AGE_FULL,3)]);
lon_vals = repmat(lon_vals,[1, size(AGE_FULL,2), size(AGE_FULL,3)]);

lon_flat = lon_vals(:);
lat_flat = lat_vals(:);
age_flat = AGE_FULL(:);

valid_indices = (lat_flat >= -35) & (age_flat >= 20);
%valid_indices = age_flat >= 20;

% Extract valid positions
valid_lon = lon_flat(valid_indices);
valid_lat = lat_flat(valid_indices);

% Top

croco_lon_flat = croco_lon(:);
croco_lat_flat = croco_lat(:);
croco_top_flat = croco_top(:);

% Step 2: Find nearest depths
valid_depths = nan(size(valid_lon)); % Preallocate for efficiency

for i = 1:length(valid_lon)
    % Find the nearest grid point for each valid particle
     distances = sqrt((croco_lon_flat - valid_lon(i)).^2 + (croco_lat_flat - valid_lat(i)).^2);
     % Find the index of the minimum distance
     [~, min_idx] = min(distances);
    % Extract corresponding depth
    valid_depths(i) = croco_top_flat(min_idx);
end

% Step 3: Compute the mean depth (ignoring NaNs)
mean_depth = mean(valid_depths, 'omitnan');
median_depth = median(valid_depths, 'omitnan');
std_depth = std(valid_depths, 'omitnan');

disp(['Mean depth of particles with age >= 20 and north of -35S: ', num2str(mean_depth), ' m']);

%% Question 2

% What is the percentage of floats greater than 20 days old initiated north of 35 S and how many of them are there  

% Step 1: Remove NaN values from all arrays
valid_idx = ~isnan(lon_flat) & ~isnan(lat_flat) & ~isnan(age_flat);

lon_flat = lon_flat(valid_idx);
lat_flat = lat_flat(valid_idx);
age_flat = age_flat(valid_idx);

% Step 2: Apply conditions
valid_points = (lat_flat >= -35) & (age_flat >= 20);


% Step 3: Calculate percentage
percentage_valid = (sum(valid_points) / length(lon_flat)) * 100;

% Display results
disp(['Percentage of points meeting conditions: ', num2str(percentage_valid), '%']);

%% Question 3

% Percentage of floats north of 35 S that are 20 days or older relative to
% total released

pn = double(ncread(strcat(CROCO_path,'/',CROCO_file),'pn'));
pm = double(ncread(strcat(CROCO_path,'/',CROCO_file),'pm'));

pn = mean(mean(pn));
pm = mean(mean(pm));

grid_area = 1/(pn*pm);  % m^2

for i = 1:size(AGE_FULL,3)
    for j = 1:size(AGE_FULL,2)
        tmp = AGE_FULL(:,j,i);
        mlat = lat_vals(:,j,i);
        mlon = lon_vals(:,j,i);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        valid_idx = ~isnan(mlon) & ~isnan(mlat) & ~isnan(tmp);
        mlat = mlat(valid_idx);
        mlon = mlon(valid_idx);
        tmp = tmp(valid_idx);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 2: Apply conditions
        valid_points = (mlat >= -35) & (tmp >= 20);
        store(j,i) = sum(valid_points);
    end
end

% Calculate percentage
num_35_full = sum(float_lat(:,1,1) >= -35);  % Number of initiated floats

out_per = mean(store(:)./num_35_full);
out_std = std(store(:)./num_35_full);

% Display results
disp(['Percentage of points meeting conditions relative to total: ', num2str(out_per), '+-',num2str(out_std)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute MLD flux for the region

flux = store.*grid_area;  % m^2
flux_mean = mean(flux(:));
flux_std_mean = std(flux(:));

% Equation of fit is xxxx for export

ff_fit = mean(mean_NSINK(20:30));  % mean value of sinking flux between 20-30 days in units of mmol C m^-2 day^-1

my_flux_days = (flux_mean * ff_fit)*0.01201; % Convert to units of gC/day
my_flux_std = (flux_std_mean * ff_fit)*0.01201;

% Compute the export flux
mean_NEXP = mean_NEXP * CtoN; % Convert units

ff_fit = mean(mean_NEXP(20:30));  % mean value of export flux between 20-30 days in units of mmol C m^-2 day^-1

my_exp_days = (flux_mean * ff_fit)*0.01201; % Convert to units of gC/day
my_exp_std = (flux_std_mean * ff_fit)*0.01201;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[myVAR, age_vals, lon_vals, lat_vals]= agetovar(float_lon, float_lat, age_matrix(:,ind_float),5, CHL(:,:,ind_chl), croco_lon, croco_lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPARE

% Chl and Sinking flux: decay function
% figure 
% plot(1:30,mean_CHL,'o')
% curve = showfit('f(t)=a*exp(lambda*t)+b','fitlinewidth',2,'fitcolor','red');
%         title(strcat('Best fit for data'));
% 
% % DET: Polynomial function
% figure 
% plot(1:30,mean_NSINK.*-1,'o')
% curve = showfit('f(t)= poly2','fitlinewidth',2,'fitcolor','red');
%         title(strcat('Best fit for data'));
