%%%% Main data extraction script for CROCO fluxes
%%%% This script takes a flux input variable from a CROCO file and converts
%%%% it to units of mmol N/m^2/s for a specific depth chosen by the user
%%%% In brief, the script converts the volume flux to an area flux by
%%%% multiplying the variable by the layer depth.
%%%% Only interested in January. The correspinding period of the Expriments
%%%%
%%%% Flux variables that can be extracted include:
%%%% New production - NFlux_NewProd - NPROD
%%%% Zooplankton grazing - NFlux_Grazing - NGRAZ
%%%% Zooplnakton mortality -  NFlux_Zmort - NZOO
%%%% Remineralisation rate - NFlux_ReminD - REM
%%%% POC sinking flux - NFlux_VSinkD1 - NSINK
%%%% Carbon export (CMIP) - NFlux_VSinkD1 - NEXP - depth horizon is 100 m
%%%% for continuity with other scripts, use suggested variable names
%%%% units are mmol N m-3 s-1
%%%% Important to note that NSINK is on the w grid and the other fluxes on
%%%% the rho, this needs to be taken into account
%%%%
%%%% Written: Jonathan Rogerson
%%%% June 2025
%%%% NOTE: CORRECTED version. To avoid fuzzy logic, interpolation is first
%%%% done to depth level and then the layer thickness is multiplied.

addpath('/usr/local/MATLAB/R2025b/toolbox/matlab/imagesci_utils/')

CROCO_path = '/media/jrogerson/JONO/CROCO_BIO';
Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

% SEASON
% Summer is January and winter is July. 
MS = 7;      % 1 for Summer and 7 for winter

% User inputs
var = 'NFlux_VSinkD1'; % variable to process
var_type = 'NSINK'; % output name
type = 'w'; % rho grid or w-grid
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
CROCO_top = ncread(file,'h');

%%%%%%%%%%%%%%%%%%%%%%%%%
% DO an integrate
%%%%%%%%%%%%%%%%%%%%%%%%%

start
% 
%  !!! WARNING weak point: vtransform should be the one used for CROCO
%
vtransform=2;
%
% Levele
level = -30;  % Depth to extract (default is -30)

%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%

VAR_mean = [];
TIME = [];

for yr = Ymin:Ymax

    %--------------------------------------------------------------
    % Determine which months to read
    % Summer : December (previous year) + January (current year)
    % Winter : June + July (current year)
    %--------------------------------------------------------------

    if MS == 1
        years_to_read  = [yr-1, yr];
        months_to_read = [12, 1];
    elseif MS == 7
        years_to_read  = [yr, yr];
        months_to_read = [6, 7];
    else
        error('MS must be 1 (summer) or 7 (winter).')
    end

    %--------------------------------------------------------------
    % Loop over the required months
    %--------------------------------------------------------------

    for jj = 1:2

        % CROCO diabio file
        fname = strcat(CROCO_path,'/','diabio_avg_Y', ...
            string(years_to_read(jj)),'M',string(months_to_read(jj)),'.nc');
        fname = convertStringsToChars(fname);

        disp(fname)
        disp('PROCESSING')

        % Read time vector
        nc = netcdf(fname);
        time = nc{'time'}(:);
        close(nc)

        %----------------------------------------------------------
        % Loop over all timesteps
        %----------------------------------------------------------

        for l = 1:length(time)

            % Read diabio variables
            nc = netcdf(fname);

            h       = nc{'h'}(:);
            theta_s = nc.theta_s(:);
            theta_b = nc.theta_b(:);
            hc      = nc.hc(:);
            N       = length(nc('s_rho'));

            myvar = squeeze(nc{var}(l,:,:,:));

            % Convert w-grid variables to rho-grid
            if type == 'w'
                myvar = 0.5 * (myvar(1:end-1,:,:) + myvar(2:end,:,:));
            end

            close(nc);

            %------------------------------------------------------
            % Read zeta and time from the corresponding avg file
            %------------------------------------------------------

            fname_z = strrep(fname,'diabio_','');

            nc = netcdf(fname_z);

            zeta = squeeze(nc{'zeta'}(l,:,:));
            time_now = nc{'time'}(l);

            TIME = cat(1,TIME,time_now);

            close(nc);

            %------------------------------------------------------
            % Interpolate flux to chosen depth
            %------------------------------------------------------

            z = get_depths(fname,fname,l,'r');

            myvar_interp = vinterp(myvar,z,level);

            %------------------------------------------------------
            % Compute layer thickness
            %------------------------------------------------------

            zw = zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);

            dz = zw(2:end,:,:) - zw(1:end-1,:,:);

            [~,ny,nx] = size(dz);

            dz_d = nan(ny,nx);

            for m = 1:ny
                for ii = 1:nx

                    z_profile = squeeze(z(:,m,ii));

                    [~,kd] = min(abs(z_profile + abs(level)));

                    dz_d(m,ii) = dz(kd,m,ii);

                end
            end

            %------------------------------------------------------
            % Convert volumetric flux to vertically integrated flux
            %------------------------------------------------------

            out = myvar_interp .* dz_d;      % mmol N m^-2 s^-1

            VAR_mean = cat(3,VAR_mean,out);

            clear out myvar myvar_interp dz_d z zw dz

        end

        disp('Completed month')

    end

    disp(['Completed trajectory year ',num2str(yr)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to save VAR_mean for our desired variable according to the names
% specified in the header of this script so that the files will be 
% compatible with other scripts

% Format time

% [~, ia, ~] = unique(TIME);
% VAR_mean = VAR_mean(:,:,ia);
% TIME = TIME(ia);

% CHECK

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save yearly trajectory fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = Ymin:Ymax;

for y = years

    %--------------------------------------------------------------
    % Select the months corresponding to the trajectory year
    %--------------------------------------------------------------

    if MS == 1
        % Summer: December (previous year) + January (current year)
        idx = (year(CROCO_time) == y-1 & month(CROCO_time) == 12) | ...
              (year(CROCO_time) == y   & month(CROCO_time) == 1);

    elseif MS == 7
        % Winter: June + July (current year)
        idx = (year(CROCO_time) == y & month(CROCO_time) == 6) | ...
              (year(CROCO_time) == y & month(CROCO_time) == 7);

    else
        error('MS must be 1 (summer) or 7 (winter).')
    end

    %--------------------------------------------------------------
    % Extract the data for this trajectory year
    %--------------------------------------------------------------

    VAR = VAR_mean(:,:,idx);
    time_year = CROCO_time(idx);

    %--------------------------------------------------------------
    % Save the data
    %--------------------------------------------------------------

    filename = strcat(var_type,'_',string(y),'.mat');
    save(filename,'VAR','time_year','-v7.3');

    disp(['Saved ',filename,' (', ...
          datestr(time_year(1),'dd-mmm-yyyy'),' to ', ...
          datestr(time_year(end),'dd-mmm-yyyy'),')'])

end