%%%% This script is mean to take the summer chlorophyll surface field of
%%%% the mdodel and compare it to the observations from the OC-CCI.
%%%% Will first process the model data and then the OBS and get them to
%%%% overlap. Will also need to compute some metrics. Mean, standard
%%%% devation and bias (both in the near shore and 'offshore environment
%%%% Written: Jonathan Rogerson

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/Transport
addpath /usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/

CROCO_path = '/media/jono/JONO/CROCO_BIO';
Ymin = 2011;
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
CROCO_lat=double(ncread(strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc'),'lat_rho'));
CROCO_lon=double(ncread(strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(1),'.nc'),'lon_rho'));
CROCO_top = ncread(file,'h');
CROCO_top = double(CROCO_top .*mask);

CHL = [];
TIME = [];

for i = Ymin:Ymax
    for j = 1:12
        % CROCO average name
        fname=strcat(CROCO_path,'/','avg_Y',string(i),'M',string(j),'.nc');
        disp('Processing')
        disp(fname)
        % Read in the surface chlorophyll and time
        chl = double(ncread(fname,'CHLA'));
        chl = squeeze(chl(:,:,50,:)); % Only the surface layer please
        time = ncread(fname,'time');
        % SAVE the outputs
        disp('Saving CHL and time')
        CHL = cat(3,CHL,chl);
        TIME = cat(1,TIME,time);
        clear chl time   % Just for space and to be tidy
    end
end

% Process time data and mask variables as needed

% Convert the time data

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);  

%% Process the observations 

% File name to the OC-CCI file
% This file was processed using the following script and file: 
% cdo remapbil,model_grid.txt CCI_CHL_v4.0_daily_2001-2018_sub.nc obs_out.nc
% model_grid.txt was created from the CROCO_file
file = '/home/jono/Documents/POSTDOC/DATA/Roff_OUT/CCI_CHL_v4.0_daily_2001-2018_sub_INTERP.nc';

% obs_lon = double(ncread(file,'lon'));
% obs_lat = double(ncread(file,'lat'));
obs_chl = ncread(file,'chlor_a');
obs_time = ncread(file,'time');

OBS_time = datetime(1970,1,1) + days(obs_time);

% Subset

% Extract the start and end dates of CROCO_time
start_date = min(CROCO_time);  % Start of CROCO time period
end_date = max(CROCO_time);    % End of CROCO time period

% Find indices where OBS_time falls within CROCO_time range
obs_indices = find(OBS_time >= start_date & OBS_time <= end_date);
obs_indices = [obs_indices(1)-1; obs_indices];

% Subset OBS_time using the indices
OBS_time_subset = OBS_time(obs_indices);

% Subset obs_chl to the same dimensions
CHL_OBS = obs_chl(:,:,obs_indices);

% Average all variable data for our seasons
summer = [1, 2, 12]; 

ind_sum = [];
for i = 1:length(summer)
    ind_sum = [ind_sum; find(summer(i) == MO)];
end

CROCO_CHL_sum = mean(CHL(:,:,ind_sum),3,'omitnan');   % Model
OBS_CHL_sum = mean(CHL_OBS(:,:,ind_sum),3,'omitnan');

% Compute RMSE

RMSE = sqrt(mean((CHL - CHL_OBS).^2,3,'omitnan')); % Root Mean Squared Error
BIAS = CROCO_CHL_sum - OBS_CHL_sum;

% Normalize
% OBS_data = OBS_CHL_sum./(max(OBS_CHL_sum(:)-min(OBS_CHL_sum(:))));
% MODEL_data = CROCO_CHL_sum./(max(CROCO_CHL_sum(:)-min(CROCO_CHL_sum(:))));

%% Figure design:
% The figure will be designed as such with 4 subplots in a 2x2 fashion
% Subplot 1: Model CHL summer field with small map inlet of Southern Africa
% Subplot 2: OBS CHL summer filed (subplot 1 and 2 will share the same
% colorbar
% Subplot 3: Bias plot with the balance colorbar: spatial mean and std for
% model and obs must be shown in a small legend
% Subplot 4: RMSE plot with legend showing mean RMSE value for domain

myposition  = [675,123,684,832];

% Range for CHL
cmin = 0;
cmax = 10;
levels = [cmin:0.5:cmax];
mydepths = [100,200,300];

f = figure;
f.Position = myposition;
subplot(2,2,1) % MODEL CHLORPHYLL field
m_proj('miller','long',[min(min(CROCO_lon)), max(max(CROCO_lon))], ...
            'lat',[min(min(CROCO_lat)), max(max(CROCO_lat))]);
m_contourf(CROCO_lon,CROCO_lat,CROCO_CHL_sum,levels);
shading flat
cmocean('algae',length(levels))
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
    [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
title('(a) CROCO [DJF]','fontsize',15);
caxis([cmin cmax])
cRange=caxis;
hold on
[C,h] = m_contour(CROCO_lon,CROCO_lat,CROCO_top,[200 200],'ShowText','on');
h.LineWidth = 3;
h.Color = 'r';
caxis(cRange)
% In lay
axes('position',[0.32 0.81 0.12 0.12])
box on
% Define the region of interest
lon_min = min(CROCO_lon(:)); lon_max = max(CROCO_lon(:));
lat_min = min(CROCO_lat(:)); lat_max = max(CROCO_lat(:));
% Create a basic map using m_map
m_proj('mercator', 'longitudes', [10 40], 'latitudes', [-40 -15]); % Zoomed-out view of southern Africa
m_coast('line', 'Color', 'k'); % Add coastlines
m_grid('linewi', 1, 'xticklabels',[],'yticklabels',[],'tickdir', 'out', 'fontsize', 10, 'box', 'fancy'); % Add grid
% Draw the region of interest (box)
hold on;
m_line([lon_min lon_max lon_max lon_min lon_min], ...
       [lat_min lat_min lat_max lat_max lat_min], ...
       'color', 'k', 'linewidth', 2);
m_gshhs_i('line', 'Color', 'k', 'LineStyle', '--');

subplot(2,2,2) % OBS CHL field
m_proj('miller','long',[min(min(CROCO_lon)), max(max(CROCO_lon))], ...
            'lat',[min(min(CROCO_lat)), max(max(CROCO_lat))]);
m_contourf(CROCO_lon,CROCO_lat,OBS_CHL_sum,levels);
shading flat
cmocean('algae',length(levels))
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
    [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
title('(b) OBS [DJF]','fontsize',15);
caxis([cmin cmax])
cRange=caxis;
hold on
[C,h] = m_contour(CROCO_lon,CROCO_lat,CROCO_top,[200 200],'ShowText','on');
h.LineWidth = 3;
h.Color = 'r';
caxis(cRange)
% Create colorbar
hold on
ca = colorbar('southOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = '[mg Chla m-3]';
ca.FontSize = 11;
ca.Position = [0.377360731602555 0.533642691415313 0.301313619982458 0.0108584693011625]; 
caxis([cmin cmax]);

subplot(2,2,3) % BIAS PLOT
cmin = -5;
cmax = 5;
levels = [cmin:0.5:cmax];
m_proj('miller','long',[min(min(CROCO_lon)), max(max(CROCO_lon))], ...
            'lat',[min(min(CROCO_lat)), max(max(CROCO_lat))]);
m_pcolor(CROCO_lon,CROCO_lat,BIAS);
shading flat
cmocean('balance',length(levels))
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
    [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
title('(c) Bias [model-obs]','fontsize',15);
caxis([cmin cmax])
cRange=caxis;
hold on
[C,h] = m_contour(CROCO_lon,CROCO_lat,CROCO_top,[200 200],'ShowText','on');
h.LineWidth = 3;
h.Color = 'r';
caxis(cRange)
%
hold on
ca = colorbar('southOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = '[mg Chla m-3]';
ca.FontSize = 11;
ca.Position = [0.153081937770145 0.0603248259860789 0.286399330241382 0.0120659700245637]; 
caxis([cmin cmax]);
%
% Text-box
% Add text box with border and background color
% Caculate mean and std for the model and obs
mean_model = mean(CROCO_CHL_sum(:),'omitnan');
std_model = std(CROCO_CHL_sum(:),'omitnan');
mean_obs = mean(OBS_CHL_sum(:),'omitnan');
std_obs = std(OBS_CHL_sum(:),'omitnan');
hold on
stat_text = {'mean ± std', 'CROCO: 0.89 ± 1.11','OBS: 1.90 ± 2.44'};
annotation('textbox', [0.305 0.345 0.1 0.1], ... % Position: [x y width height] (normalized)
           'String', stat_text, ...
           'FontSize', 8, ...
           'BackgroundColor', [.7 .7 .7], ... % White background
           'EdgeColor', 'k', ... % Black border
           'LineWidth', 1.5, ...
           'HorizontalAlignment', 'center');

subplot(2,2,4) % This is the RMSE
cmin = 0;
cmax = 10;
levels = [cmin:0.5:cmax];
m_proj('miller','long',[min(min(CROCO_lon)), max(max(CROCO_lon))], ...
            'lat',[min(min(CROCO_lat)), max(max(CROCO_lat))]);
m_pcolor(CROCO_lon,CROCO_lat,RMSE);
shading flat
cmocean('matter',length(levels))
hold on
m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
    [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
title('(d) RMSE','fontsize',15);
caxis([cmin cmax])
cRange=caxis;
hold on
[C,h] = m_contour(CROCO_lon,CROCO_lat,CROCO_top,[200 200],'ShowText','on');
h.LineWidth = 3;
h.Color = 'r';
caxis(cRange)
%
hold on
ca = colorbar('southOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = '[mg Chla m-3]';
ca.FontSize = 11;
ca.Position = [0.593499060132815 0.0614849187935035 0.286904398080442 0.0108584693011631]; 
caxis([cmin cmax]);
%
% Text-box
% Add text box with border and background color
hold on
stat_text = {'Mean RMSE: 2.39'};
annotation('textbox', [0.73 0.34 0.1 0.1], ... % Position: [x y width height] (normalized)
           'String', stat_text, ...
           'FontSize', 8, ...
           'BackgroundColor', [.7 .7 .7], ... % White background
           'EdgeColor', 'k', ... % Black border
           'LineWidth', 1.5, ...
           'HorizontalAlignment', 'center');

% Save the figure
set(gcf, 'InvertHardcopy', 'off')
print('-f1','SBUS_Fig1','-dpng','-r600');






































