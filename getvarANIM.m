%%%% This script will create animations of a chosen tracer for model outputs 
% Written: Jonathan Rogerson
% Date: June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIVERSITY of CAPE TOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script calls CROCO outputs 
% Objective is to animate over variable of choice fo the model period
% Remember, visual outputs are very expensive 

%% Get the main File paths for the input data

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/Transport
addpath /usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/

CROCO_path = '/media/jono/JONO/CROCO_BIO';
Ymin = 2011;
Ymax = 2011;
Yorig = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','avg_Y',string(Ymin),'M',string(6),'.nc');
% Read in the mask
mask=ncread(file,'mask_rho');
mask(mask==0)=nan;
% Define my region, will focus on the coastal domain
CROCO_lat=ncread(strcat(CROCO_path,'/','diabio_avg_Y',string(Ymin),'M',string(6),'.nc'),'lat_rho');
CROCO_lon=ncread(strcat(CROCO_path,'/','diabio_avg_Y',string(Ymin),'M',string(6),'.nc'),'lon_rho');
CROCO_top = ncread(file,'h');

var = ncread(file,'CHLA');
var = squeeze(var(:,:,50,:));

time = ncread(file,'time');
CROCO_time = datetime(Yorig,1,1) + seconds(time);
[Y,MO,D] = datevec(CROCO_time); 

% Declare double

CROCO_lon = double(CROCO_lon);
CROCO_lat = double(CROCO_lat);
CROCO_top = double(CROCO_top);
var = double(var);

%% Figure and animation

cmin = 0;
cmax = 5;
levels = [cmin:0.5:cmax];
mydepths = [100,200,300];

for i = 1:length(time)
    clf
    figure(1); 
    filename = 'CROCO_anim_sum.gif';
    hold on;
    title(string(CROCO_time(i)),'fontsize',15);

    % Plot loop
    m_proj('miller','long',[min(min(CROCO_lon)), max(max(CROCO_lon))], ...
                'lat',[min(min(CROCO_lat)), max(max(CROCO_lat))]);
    m_pcolor(CROCO_lon,CROCO_lat,var(:,:,i));
    shading flat
    cmocean('algae',length(levels))
    hold on;
    m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
        [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
    m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    caxis([cmin cmax])
    cRange=caxis;
    ca = colorbar('southOutside');
    ca.Label.String = '[mg Chla m-3]';
    ca.FontSize = 11;
    caxis([cmin cmax]);
    % Animate and save
    pause(0.02);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
         imwrite(imind,cm,filename,'gif','WriteMode','append');
    end 
end








