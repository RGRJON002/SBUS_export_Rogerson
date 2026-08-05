%%%% This script takes the summer (january) and winter (july) release
%%%% experiments and uses the lon, lat and time to index the eulerian
%%%% fields extracted used the two NPZD_get functions. The goal is to save
%%%% a sttructure for summer and winter, containing the age and assoicated
%%%% mean tracer evolution + std. For this script, the goal is to look at 5
%%%% tracers: NO3, CHLA, DET, NSINK and REM.
%%%%
%%%% Written: Jonathan Rogerson
%%%% This code builds off an earlier version, but fixes the sampling regime
%%%% for the file

addpath('/usr/local/MATLAB/R2025b/toolbox/matlab/imagesci_utils/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File paths to data
float_path = '/home/jrogerson/Documents/AJMS_project/DATA/Roff_OUT';
CROCO_path = '/media/jrogerson/JONO/CROCO_BIO/'; 
CROCO_file = strcat(CROCO_path,'avg_Y2012M1.nc');

Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

Fcount = 31;

% NB

% My age tracer to exclude

cut_age = 0;

% SEASON
% Summer is January and winter is July. 
MS = 7;      % 1 for Summer and 7 for winter

% MAT files: The location of my summer and winter .mat fields
if MS == 1 % SUMMER
    mat_path = '/home/jrogerson/Documents/AJMS_project/MANUSCRIPT/Version_revision/Figures/Figure_lagrangain fluxes/SUMMER/';
elseif MS == 7   % WINTER
    mat_path = '/home/jrogerson/Documents/AJMS_project/MANUSCRIPT/Version_revision/Figures/Figure_lagrangain fluxes/WINTER/';
else
    disp('NOT A VALID ENTRY')
end

% Want to save and store in a stacked array the following processed
% variables
RESULTS = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = Ymin:Ymax
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: Process the Lagrangian files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('STEP 1: Processing the Lagrangian files')

    % DECLARE the files to loop over
    %float file
    if MS == 1
        disp('PROCESSING SUMMER')
        float_file = strcat('SBUS_floats_bio_Jan',string(i),'.nc');
        disp(float_file)
    elseif MS == 7
        disp('PROCESSING WINTER')
        float_file = strcat('SBUS_floats_bio_Jul',string(i),'.nc');
        disp(float_file)
    else
    end

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
    % Convert 6-hourly timestamps to daily timestamps
    time_float = FUtime(time_float);

    % Clean the files
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % avgerage over days
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    float_lon = FU3D(float_lon);
    float_lat= FU3D(float_lat);
    float_depth = FU3D(float_depth);
    float_temp = FU3D(float_temp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECTION 2: READ in the .mat files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % file names of .mat
    
    % NSINK file
    SINK_file = strcat(mat_path,'NSINK_',string(i),'.mat');
    % DET file
    DET_file = strcat(mat_path,'DET_',string(i),'.mat');
    % REM file 
    REM_file = strcat(mat_path,'REM_',string(i),'.mat');
    % CHLA
    CHLA_file = strcat(mat_path,'CHLA_',string(i),'.mat');
    % NO3
    NO3_file = strcat(mat_path,'NO3_',string(i),'.mat');

    %%%%%%%%%%%%
    % Load CROCO
    %%%%%%%%%%%%
    
    disp('READING CROCO grid information')
    % Read in the grid dimension of one of the CROCO_files
    
    croco_lon = double(ncread(CROCO_file,'lon_rho'));
    croco_lat = double(ncread(CROCO_file,'lat_rho'));
    
    % Read in the topography data

    croco_top = ncread(CROCO_file,'h');
    mask = ncread(CROCO_file,'mask_rho');
    mask(mask==0)=nan;               % Land is = 0, make nan
    croco_top=croco_top.*mask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('PROCESSING .mat files')
    % Read in the .mat files. We can access it using dot(.) indexing
     % NSINK
    nsink = load(SINK_file);
    nsink.VAR = permute(nsink.VAR, [2, 1, 3]);  % Get the dimensions aligned with croco_lon,lat ect...
    nsink.VAR = nsink.VAR.*mask;        % Mask as you would any variable
    %%%%%%%%%%%%%%%%%%%%%
    % Once off time read
    %%%%%%%%%%%%%%%%%%%%%
    time_croco = nsink.time_year;

    % DET
    det = load(DET_file);
    det.VAR = permute(det.VAR, [2, 1, 3]);
    det.VAR = det.VAR.*mask;
    
    % NREM
    rem = load(REM_file);
    rem.VAR = permute(rem.VAR, [2, 1, 3]);
    rem.VAR = rem.VAR.*mask;

    % CHLA
    chla = load(CHLA_file);
    chla.VAR = permute(chla.VAR, [2, 1, 3]);
    chla.VAR = chla.VAR.*mask;

    % NO3
    no3 = load(NO3_file);
    no3.VAR = permute(no3.VAR, [2, 1, 3]);
    no3.VAR = no3.VAR.*mask;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STAGE 3: Lagrangian sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % We also want to sample bathymetry so will do a small hack and repmat
    % the topography to the third dimension size of the croco_time :)
    
     % Depth
    dpt = repmat(croco_top, [1 1 size(time_croco,1)]);

    disp('LAGRNGAIN SAMPLING')

    NSINK = sample_croco( ...
        nsink.VAR, croco_lon, croco_lat, time_croco, ...
        float_lon, float_lat, time_float);

    DET = sample_croco( ...
        det.VAR, croco_lon, croco_lat, time_croco, ...
        float_lon, float_lat, time_float);
    
    REM = sample_croco( ...
       rem.VAR, croco_lon, croco_lat, time_croco, ...
       float_lon, float_lat, time_float);
    
    CHLA = sample_croco( ...
        chla.VAR, croco_lon, croco_lat, time_croco, ...
        float_lon, float_lat, time_float);
    
    NO3 = sample_croco( ...
        no3.VAR, croco_lon, croco_lat, time_croco, ...
        float_lon, float_lat, time_float);

    DPT = sample_croco( ...
        dpt, croco_lon, croco_lat, time_croco, ...
        float_lon, float_lat, time_float);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 4: Age MAtrix computation and masking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Our sampled trajectories at the moment include just 30 days backwards
    % However, we want to track/truncate the temporal tracking until
    % particles reach the 200 m isobath, which we use as our coastal
    % domain/origin. We can leverage the computation of the age_matrix and
    % use it as a mask for our other variables: IE - We track the evolution
    % till partilces reach the 200 m isobath

    disp('AGE MATRIX')

    coastal_age = AGE_TRACER(DPT);

% ==========================================================
% DIAGNOSTICS
% ==========================================================

    fprintf('\n');
    fprintf('Particles reaching the coast\n');
    fprintf('----------------------------\n');
    
    fprintf('Valid trajectories : %d of %d\n', ...
        sum(~isnan(coastal_age(:))), numel(coastal_age));
    
    fprintf('Fraction reaching coast : %.1f %%\n', ...
        100*sum(~isnan(coastal_age(:)))/numel(coastal_age));
    
    fprintf('Mean coastal age : %.2f days\n', ...
        nanmean(coastal_age(:)));
    
    fprintf('Median coastal age : %.2f days\n\n', ...
        nanmedian(coastal_age(:)));

    NSINK = mask_by_age(NSINK,coastal_age);
    DET   = mask_by_age(DET,coastal_age);
    REM   = mask_by_age(REM,coastal_age);
    CHLA  = mask_by_age(CHLA,coastal_age);
    NO3   = mask_by_age(NO3,coastal_age);
    DPT   = mask_by_age(DPT,coastal_age);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 5: Compute the age-based statisitics :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('COMPUTING AGE-BASED STATISTICS')

    NSINK_stats = age_statistics(NSINK);
    DET_stats   = age_statistics(DET);
    REM_stats   = age_statistics(REM);
    CHLA_stats  = age_statistics(CHLA);
    NO3_stats   = age_statistics(NO3);
    DPT_stats   = age_statistics(DPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE 6: Save the statistics for each :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('SAVING')

    RESULTS(i-Ymin+1).year = i;

    RESULTS(i-Ymin+1).NSINK = NSINK_stats;
    RESULTS(i-Ymin+1).DET   = DET_stats;
    RESULTS(i-Ymin+1).REM   = REM_stats;
    RESULTS(i-Ymin+1).CHLA  = CHLA_stats;
    RESULTS(i-Ymin+1).NO3   = NO3_stats;
    RESULTS(i-Ymin+1).DPT   = DPT_stats;


   disp(['FINISHED with ', string(i)])

   clear NSINK_stats DET_stats CHLA_stsats NO3_stats DPT_stats REM_stats
            
end

disp('FINISHED')

%% SAVE

% We will process the outputs in a .py script for the plotting 

save('Winter_AgeStatistics.mat','RESULTS','-v7.3')







