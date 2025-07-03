%%%% This script is intended to visualize the initial position of the
%%%% floats along with the seasonal temperture and geostrophic currents

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/Documents/PhD/CHAPTER_2/DATA
addpath /home/jono/Documents/PhD/CHAPTER_2/SCRIPTS

%float file
float_file = 'SBUS_floats_bio_Jan2012.nc';
float_path = '/home/jono/Documents/POSTDOC/DATA/Roff_OUT';

% Inital file
CROCO_path = '/media/data/Roff'; 
CROCO_file = 'croco_N01.nc';

% How many floats per release point

Fcount = 31;

% Depths I want to process 

Depth1 = [-5];  % My subsurface depths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Small additional process

float_lat = squeeze(float_lat(:,1,1));
float_lon = squeeze(float_lon(:,1,1));

% Read in the grid dimension of one of the CROCO_files

croco_lon = double(ncread(strcat(CROCO_path,'/',CROCO_file),'lon_rho'));
croco_lat = double(ncread(strcat(CROCO_path,'/',CROCO_file),'lat_rho'));

% Read in the topography data

croco_top = ncread(strcat(CROCO_path,'/',CROCO_file),'h');
mask = ncread(strcat(CROCO_path,'/',CROCO_file),'mask_rho');
mask(mask==0)=nan;               % Land is = 0, make nan
croco_top=croco_top.*mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot initial distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some variables

mydepths = [100,200,300,500];

figure
m_proj('miller','long',[min(min(croco_lon)), max(max(croco_lon))], ...
        'lat',[min(min(croco_lat)), max(max(croco_lat))]);
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
hold on
m_contour(croco_lon,croco_lat,croco_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
    hold on
    for i = 1:size(float_lat,1)
        m_plot(float_lon(i,1),float_lat(i,1), 'ko','MarkerSize',3,'MarkerFaceColor','g') %plot the floats
        hold on
    end
hold on
[C,h] = m_contour(croco_lon,croco_lat,croco_top,[300 300],'ShowText','on');
h.LineWidth = 3;
h.Color = 'r';
hold on

% Save the figure

set(gcf, 'InvertHardcopy', 'off')
print('-f1','init_figure_bii','-dpng','-r600');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







