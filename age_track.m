%%%% Algorithm to calculate age tracer of particles. Seeded particles are
%%%% tracked backwards in time from their release location. The objective
%%%% is to see whether any particles have their origin in the coastal
%%%% waters and are exported offshore. 
%%%% To accomplish this, we are going to make use of the fact we will only
%%%% consider particles having a coastal origin if along their path they
%%%% touch or become shallower than 200 m. Else, they will be given a value
%%%% of NAN

%%%% Written: Jonathan Rogerson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/Documents/PhD/CHAPTER_2/DATA
addpath /home/jono/Documents/PhD/CHAPTER_2/SCRIPTS

%float file
float_path = '/home/jono/Documents/POSTDOC/DATA/Roff_OUT';

% Inital file
CROCO_path = '/media/jono/JONO/CROCO_BIO'; 

Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

% How many floats per release point

Fcount = 31;
cut_age = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Want to save and store in a stacked array the following processed
% variables
myCHL = [];    % 3-D
DET = [];    % 3-D
AGE_FULL=[]; % 3-D
TIME = [];   % 1-D

for i = Ymin:Ymax
    % DECLARE the files to loop over
    %float file
    float_file = strcat('SBUS_floats_bio_Jan',string(i),'.nc');
    % CROCO file
    CROCO_file = strcat('avg_Y',string(i),'M1.nc');
    % DET file
    DET_file = strcat('DET_',string(i),'.mat');
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
    TIME = cat(1,TIME,time);   % Save to array
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
    myCHL = cat(3,myCHL,CHL);
   
    % DET
    load(strcat('DET_',string(i),'.mat'));
    det = VAR;
    det = permute(det, [2, 1, 3]);
    det = det.*mask;
    DET = cat(3,DET,det);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the age
    
    age_matrix = AGE_TRACER(float_lon,float_lat, croco_lon, croco_lat, croco_top);
    
    AGE_FULL = cat(3,AGE_FULL,age_matrix);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(strcat('DONE with year --', string(i)));
    
end

%% POST Processing

% CROCO time
CROCO_time = datetime(2000,1,1) + seconds(TIME);

% Reshape arrays for structuring
CROCO_time = reshape(CROCO_time,31,4);
myCHL = reshape(myCHL, 102, 202, 31, 4);
DET = reshape(DET, 102, 202, 31, 4);

% Convert units

DET = DET .* 6.625;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a figure

% Chose a time for floats, remeber, the time-stepping isreversed
% Create the arrays for plotting in the subplots

% My SBUS mask

SBUS_mask = croco_top >= 200;
SBUS_mask = double(SBUS_mask);
SBUS_mask(SBUS_mask == 0) = NaN;

ind_year = [1,2,3,4];
ind_float = [24,13,9,3];
ind_chl = ind_float + 31 - (2*ind_float) + 1;

% Using the inital lon_lat pairs, plot as a 3D scatter plot

mydepths = [100,200,300];
cpmin=5;   % Custom
cpmax=30;   % Custom
crange= [cpmin cpmax];
skp = 1;

figure
% Column 1: AGE tracer field
for i = 1:length(ind_year)
    subplot(3,4,i)
    blub = AGE_FULL(:,ind_float(i),ind_year(i));
    m_proj('miller','long',[min(min(croco_lon)), max(max(croco_lon))], ...
            'lat',[min(min(croco_lat)), max(max(croco_lat))]);
    m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');   
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',11);
    hold on
    m_contour(croco_lon,croco_lat,croco_top,mydepths,'Color','g','LineWidth',1,'LineStyle','--');
    hold on
    [x, y] = m_ll2xy(squeeze(float_lon(:,1,1)),squeeze(float_lat(:,1,1))); % Convert coordinates to map projection
    scatter(x(1:skp:end), y(1:skp:end), 20, blub(1:skp:end), 'filled'); % Scatter plot
    % titles
    title_str = datestr(CROCO_time(ind_chl(i),ind_year(i)), 'dd-mmm-yyyy');
    title(title_str,'fontsize',14)
    % Adjust colormap
    colormap(cmocean('matter')); % Or another suitable colormap like parula, viridis, etc.
    caxis(crange)
    hold on
    [C,h] = m_contour(croco_lon,croco_lat,croco_top,[200 200],'ShowText','on');
    h.LineWidth = 3;
    h.Color = 'r';
    caxis(crange)
    hold on
end

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = 'Age Tracer [Days]';
ca.FontSize = 12;
caxis([cpmin cpmax]);

% Title
%title('(a) Tracer Age','fontsize',15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Corresponding chlorophyll map

cmin = 0;
cmax = 3;
mydepths = [100,200,300];

% Summer
for i = 1:length(ind_year)
    subplot(3,4,4+i)
    % Plot figure
    m_proj('miller','long',[min(min(croco_lon)), max(max(croco_lon))], ...
                'lat',[min(min(croco_lat)), max(max(croco_lat))]);
    m_pcolor(croco_lon,croco_lat,...
        myCHL(:,:,ind_chl(i),ind_year(i)).*SBUS_mask);
    shading flat
    cmocean('algae')
    hold on;
    m_contour(croco_lon,croco_lat,croco_top,mydepths,'Color',...
        [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
    m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',11);
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    [C,h] = m_contour(croco_lon,croco_lat,croco_top,[200 200],'ShowText','on');
    h.LineWidth = 3;
    h.Color = 'r';
    caxis(cRange)
end
    
% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = 'Chlorophyll [mg Chla m-3]';
ca.FontSize = 12;
caxis([cmin cmax]);

%title('(b) Chlorophyll','fontsize',15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Corresponding chlorophyll map

cmin = 0;
cmax = 12;
mydepths = [100,200,300];

% Summer
for i = 1:length(ind_year)
    subplot(3,4,8+i)
    % Plot figure
    m_proj('miller','long',[min(min(croco_lon)), max(max(croco_lon))], ...
                'lat',[min(min(croco_lat)), max(max(croco_lat))]);
    m_pcolor(croco_lon,croco_lat,...
        DET(:,:,ind_chl(i),ind_year(i)).*SBUS_mask);
    shading flat
    cmocean('dense')
    hold on;
    m_contour(croco_lon,croco_lat,croco_top,mydepths,'Color',...
        [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
    m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',11);
    caxis([cmin cmax])
    cRange=caxis;
    hold on
    [C,h] = m_contour(croco_lon,croco_lat,croco_top,[200 200],'ShowText','on');
    h.LineWidth = 3;
    h.Color = 'r';
    caxis(cRange)
end

% Create colorbar
hold on
ca = colorbar('eastOutside');
ca.Position = ca.Position + 1e-10;
ca.Position(2) = ca.Position(2) - 0.12;
ca.Label.String = 'POC [mmol C m-3]';
ca.FontSize = 12;
caxis([cmin cmax]);

%title('(c) POC','fontsize',15);

% SAVE
%%

set(gcf, 'InvertHardcopy', 'off')
print('-f1','AGE_snap','-dpng','-r600');








