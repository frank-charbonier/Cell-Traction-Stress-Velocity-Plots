function compute_cellvel(domainname, DICname, cellvel_savename, pix_size, time_increment, plot_radial, thr)
% uncomment below to use default function arguments
% arguments
%     % Name of domain. This is where cells are located. Set to [] if no domain
%     domainname = 'domain.tif';
%     % Name of displacment data to load
%     DICname = 'cells_DIC_results.mat';
%     % cellvel_savename
%     cellvel_savename = 'cellvel_processed.mat';
%     % Pixel size [microns]
%     pix_size = 1.3;
%     % Time between images
%     time_increment = 10; % min
%     % Assay format
%     % Set to 0 to plot x and y components, otherwise plot radial and tangential
%     plot_radial = 0;
%     % Reject displacements larger than give threshold, units are um/min
%     thr = 1;
% end
%
% COMPUTE_CELLVEL Compute cell velocities without plotting
%
% First run digital image correlation on images of the cells. Since you are
% using image correlation, which correlates a subset, cells must be
% confluent.
%
% This script requires a mat file containing the following:
%   x, y: 2D arrays containing the gridpoints on which the DIC was
%         computed. This can be made using Matlab's meshgrid command.
%         Units: pix.
%   u, v: Incremental displacements of cell image computed by image
%         correlation in horizontal and vertical directions. These are 2D
%         or 3D arrays of size (M, N, P) where M and N are the number of
%         rows and columns, which must match the size of x and y. Variable
%         P corresponds to different time points. If there is only one time
%         point, then the array is 2D (i.e., P=1).
%
% Written by Jacob Notbohm, Univerity of Wisconsin-Madison, 2015-2020
% Adapted by Frank Charbonier, Stanford University, 2023

% clear;
close all;
% clc;

%% --- COMPUTE CELL VELOCITIES ---
% Load variables associated with cell displacements and rename so they aren't overwritten
DIC_data = load(DICname, "x", "y", "u", "v", "d0");
u_cell=DIC_data.u;
v_cell=DIC_data.v;
x_cell=DIC_data.x;
y_cell=DIC_data.y;
d0 = DIC_data.d0;

% Convert from pix to um
x_cell=x_cell*pix_size;     y_cell=y_cell*pix_size;
u_cell=u_cell*pix_size;     v_cell=v_cell*pix_size;
% Convert from displacements to velocities
u_cell = u_cell/time_increment;
v_cell = v_cell/time_increment;

% Number of correlations
K = size(u_cell,3);

% Preallocate array for radial and tangential displacements
if (plot_radial==1)
    ut = zeros(size(u_cell));
    ur = zeros(size(u_cell));
end

% Initialize centroid to be center of image (needed for parfor loop)
x=DIC_data.x;
y=DIC_data.y;
xc = mean(x(:));
yc = mean(y(:));

% Get center from center of first domain image
if ~isempty(domainname)
    domain = imread(domainname,1);
    domain = double(domain); % Convert to double precision
    domain = domain/max(domain(:)); % Set max value to 1
    domain = logical(domain); % Convert to logical
    % Downsample domain
    % x and y grid points start at w0/2 and end w0/2 before the image ends.
    % First crop off edges so that domain matches start and end points of x
    % and y.
    domain = domain(round(min(y_cell(:))/pix_size):round(max(y_cell(:))/pix_size),round(min(x_cell(:))/pix_size):round(max(x_cell(:))/pix_size));
    domain = downsample(domain,d0); % downsample number of rows
    domain = downsample(domain',d0)'; % downsample number of cols

    stats = regionprops(domain);
    centroid = stats.Centroid;
    xc_idx = round(centroid(1));    yc_idx = round(centroid(2));

    % Convert back to um
    xc = xc_idx*d0*pix_size;
    yc = yc_idx*d0*pix_size;
else
    xc = mean(x(:));
    yc = mean(y(:));
end

% Turn off warning about temporary variables in parfor loop
% warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

% Using parallelization, runs ~100x faster
parfor k=1:K
    % Get centroid from domain
    if ~isempty(domainname)
        % DIC_data = load(DICname,"d0");
        % d0 = DIC_data.d0
        domain_k = imread(domainname,k);
        domain_k = double(domain_k); % Convert to double precision
        domain_k = domain_k/max(domain_k(:)); % Set max value to 1
        domain_k = logical(domain_k); % Convert to logical
        % Downsample domain
        % x and y grid points start at w0/2 and end w0/2 before the image ends.
        % First crop off edges so that domain matches start and end points of x
        % and y.
        domain_k = domain_k(round(min(y_cell(:))/pix_size):round(max(y_cell(:))/pix_size),round(min(x_cell(:))/pix_size):round(max(x_cell(:))/pix_size));
        domain_k = downsample(domain_k,d0); % downsample number of rows
        domain_k = downsample(domain_k',d0)'; % downsample number of cols

        stats = regionprops(domain_k);
        centroid = stats.Centroid;
        xc_idx = round(centroid(1));    yc_idx = round(centroid(2));

        % Convert back to um
        xc = xc_idx*d0*pix_size;
        yc = yc_idx*d0*pix_size;
    end
    
    % Cell velocity
    u_cell_k = u_cell(:,:,k);
    v_cell_k = v_cell(:,:,k);

    % --- Remove displacements that are too large ---
    % Remove displacements greater than threshold
    idx = (abs(u_cell_k)>thr) | (abs(v_cell_k)>thr);
    u_cell_k(idx) = nan;
    v_cell_k(idx) = nan;
    u_cell_k = inpaint_nans(u_cell_k);
    v_cell_k = inpaint_nans(v_cell_k);
    
    if ~isempty(domainname)
        % Correct for drift by finding mean of velocity outside domain
        SE = strel('disk',5,0);
        domain_dilate = imdilate(domain_k,SE);
        u_cell_k = u_cell_k - mean(u_cell_k(~domain_dilate), 'omitnan');
        v_cell_k = v_cell_k - mean(v_cell_k(~domain_dilate), 'omitnan');
        
        %         % Correct for drift by subtracting off median velocity
        %         u_cell_k = u_cell_k - median(u_cell_k(domain));
        %         v_cell_k = v_cell_k - median(v_cell_k(domain));
        
        % Set values outside domain to nan
        u_cell_k(~domain_k) = nan;
        v_cell_k(~domain_k) = nan;
        % Update values for u_cell and v_cell
        u_cell(:,:,k) = u_cell_k;
        v_cell(:,:,k) = v_cell_k;
    else
        % Subtract off median of full velocity field. This may or may not
        % work for your images.
                u_cell_k = u_cell_k - median(u_cell_k(:));
                v_cell_k = v_cell_k - median(v_cell_k(:));
    end
    
    % OPTIONAL: compute radial and angular components of velocity instead of x and y
    if (plot_radial==1)
        % Radial velocity
        % r_grid = sqrt( (x-xc).^2 + (y-yc).^2 ); % Gridpoints
        theta = atan2( (y_cell-yc), (x_cell-xc) );
    
        ur_k = u_cell_k.*cos(theta) + v_cell_k.*sin(theta);    
        % Angular velocity    
        ut_k = u_cell_k.*cos(theta+pi/2) + v_cell_k.*sin(theta+pi/2);
        % Update values for ut and ur
        ut(:,:,k) = ut_k;
        ur(:,:,k) = ur_k;
    end
end

if (plot_radial==1)
    save(cellvel_savename, 'x_cell', 'y_cell', 'u_cell','v_cell', 'ut', 'ur', "xc", "yc");
else
    save(cellvel_savename, 'x_cell', 'y_cell', 'u_cell','v_cell', "xc", "yc");
end
% outputs:

end