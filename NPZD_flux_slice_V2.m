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
%%%% for continuity with other scripts, use suggested variable names
%%%% units are mmol N m-3 s-1
%%%% Important to note that NSINK is on the w grid and the other fluxes on
%%%% the rho, this needs to be taken into account
%%%%
%%%% Written: Jonathan Rogerson
%%%% June 2025
%%%% NOTE: CORRECTED version. To avoid fuzzy logic, interpolation is first
%%%% done to depth level and then the layer thickness is multiplied.

addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath('/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/')
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/Transport

CROCO_path = '/media/jono/JONO/CROCO_BIO';
Ymin = 2012;
Ymax = 2015;
Yorig = 2000;

% User inputs
var = 'NFlux_ReminD'; % variable to process
var_type = 'REM'; % output name
type = 'r'; % rho grid or w-grid
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
level = -30;  % Depth to extract

%%%%%%%%%%%%%%
% VARIABLES
%%%%%%%%%%%%%%

VAR_mean = [];
TIME = [];

for n = Ymin:Ymax
    for j = 1
        % CROCO average name
        fname=strcat(CROCO_path,'/','diabio_avg_Y',string(n),'M',string(j),'.nc');
        fname = convertStringsToChars(fname);
        disp(fname)
        disp('PROCESSING')
        nc=netcdf(fname);
        time = nc{'time'}(:);  % Used to loop, will save Time from avg file
        close(nc)
    %
    % Time index:
        for l = 1:length(time)
    % Read data
            nc=netcdf(fname);
            h=nc{'h'}(:);
            theta_s=nc.theta_s(:);
            theta_b=nc.theta_b(:);
            hc=nc.hc(:);
            N=length(nc('s_rho'));
            myvar=squeeze(nc{var}(l,:,:,:));
            % Convert to rho gid
            if type == 'w'
                myvar = 0.5 * (myvar(1:end-1,:,:) + myvar(2:end,:,:));
            end
            close(nc);
            % ZETA and time must be obtained form the avg file and not the diabio
            fname_z = strrep(fname, 'diabio_', '');
            nc=netcdf(fname_z);
            zeta=squeeze(nc{'zeta'}(l,:,:));
            time = nc{'time'}(l);
            TIME = cat(1,TIME,time);
            close(nc);

             % 1. Get the full 3D flux (mmol N / m^3 / s)
            z=get_depths(fname,fname,l,'r');
            myvar_interp = vinterp(myvar, z, level);  % interpolates to depth, mmol N/m^3/s

            % 2. Get the layer thickness at chosen
            % Use zlevs to get the w-grid (layer interfaces), then compute dz
            zw = zlevs(h, zeta, theta_s, theta_b, hc, N, 'w', vtransform);
            dz = zw(2:end,:,:) - zw(1:end-1,:,:);
            
            % 3. Layer thickness 
            [nz, ny, nx] = size(dz);  % dz is [Nz, Ny, Nx]

            % Preallocate the output
            dz_d = nan(ny, nx);

            for m = 1:ny
                for n = 1:nx
                    % Get the full depth profile at each grid point
                    z_profile = squeeze(z(:,m,n));

                    % Find index of the level closest to 100 m
                    [~, kd] = min(abs(z_profile + abs(level)));

                    % Extract dz at that level
                    dz_d(m,n) = dz(kd,m,n);
                end
            end

            % 4. Convert volumetric flux to area flux at depth
            out = myvar_interp .* dz_d;  % mmol N/m^2/s
            
            % Save the variable           
            VAR_mean = cat(3,VAR_mean,out);
            clear out
        end
    end
    disp('Completed file')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to save VAR_mean for our desired variable according to the names
% specified in the header of this script so that the files will be 
% compatible with other scripts

% Format time

[~, ia, ~] = unique(TIME);
VAR_mean = VAR_mean(:,:,ia);
TIME = TIME(ia);

% CHECK

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);

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
