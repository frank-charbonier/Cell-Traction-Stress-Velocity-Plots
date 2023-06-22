function compute_cellvel(domainname, DICname, cellvel_savename, pix_size, time_increment, plot_radial)
arguments
    % Name of domain. This is where cells are located. Set to [] if no domain
    domainname = 'domain.tif';
    % Name of displacment data to load
    DICname = 'cells_DIC_results.mat';
    % cellvel_savename
    cellvel_savename = 'cellvel_processed.mat';
    % Pixel size [microns]
    pix_size = 1.3;
    % Time between images
    time_increment = 10; % min
    % Assay format
    % Set to 0 to plot x and y components, otherwise plot radial and tangential
    plot_radial = 0;
end
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
% Load data
load(DICname);

% Rename variables associated with cell displacements so they aren't overwritten
u_cell=u; v_cell=v; x_cell=x; y_cell=y;
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

% Get center from center of first domain image
if ~isempty(domainname)
    domain1 = imread(domainname,1);
    domain1 = double(domain1); % Convert to double precision
    domain1 = domain1/max(domain1(:)); % Set max value to 1
    % Downsample domain
    % x and y grid points start at w0/2 and end w0/2 before the image ends.
    % First crop off edges so that domain matches start and end points of x
    % and y.
    domain1 = domain1(round(min(y_cell(:))/pix_size):round(max(y_cell(:))/pix_size),round(min(x_cell(:))/pix_size):round(max(x_cell(:))/pix_size));
    domain1 = downsample(domain1,d0); % downsample number of rows
    domain1 = downsample(domain1',d0)'; % downsample number of cols
    % [M, N] = size(x);
    % domain = domain(1:M,1:N); % Correct for slightly larger domain. I should clean this up later.
    % Centroid coordinates
    xc = sum(x(:).*domain1(:)) / sum(domain1(:)); % Units: pix (same units as x and y)
    yc = sum(y(:).*domain1(:)) / sum(domain1(:));
    % Find indices corresponding to nearest x and y coords to xc and yc
    xv = x(1,:);
    yv = y(:,1);
    distx = abs(xc-xv);
    disty = abs(yc-yv);
    [~, xc_idx] = min(distx);
    [~, yc_idx] = min(disty);
    % Convert back to um
    xc = xc_idx*d0*pix_size;
    yc = yc_idx*d0*pix_size;
else
    xc = mean(x(:));
    yc = mean(y(:));
end

for k=1:K
    
    % Domain
    if ~isempty(domainname)
        domain = imread(domainname,k);
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
        % Centroid coordinates
        xc = sum(x(:).*domain(:)) / sum(domain(:)); % Units: pix (same units as x and y)
        yc = sum(y(:).*domain(:)) / sum(domain(:));
        % Find indices corresponding to nearest x and y coords to xc and yc
        xv = x(1,:);
        yv = y(:,1);
        distx = abs(xc-xv);
        disty = abs(yc-yv);
        [~, xc_idx] = min(distx);
        [~, yc_idx] = min(disty);
        % Convert back to um
        xc = xc_idx*d0*pix_size;
        yc = yc_idx*d0*pix_size;
    end
    
    % Cell velocity
    u_cell_k = u_cell(:,:,k);
    v_cell_k = v_cell(:,:,k);
    
    if ~isempty(domainname)
        % Correct for drift by finding mean of velocity outside domain
        SE = strel('disk',5,0);
        domain_dilate = imdilate(domain,SE);
        u_cell_k = u_cell_k - nanmean(u_cell_k(~domain_dilate));
        v_cell_k = v_cell_k - nanmean(v_cell_k(~domain_dilate));
        
        %         % Correct for drift by subtracting off median velocity
        %         u_cell_k = u_cell_k - median(u_cell_k(domain));
        %         v_cell_k = v_cell_k - median(v_cell_k(domain));
        
        % Set values outside domain to nan
        u_cell_k(~domain) = nan;
        v_cell_k(~domain) = nan;
        % Update values for u_cell and v_cell
        u_cell(:,:,k) = u_cell_k;
        v_cell(:,:,k) = v_cell_k;
    else
        % Subtract off median of full velocity field. This may or may not
        % work for your images.
        %         u_cell_k = u_cell_k - nanmedian(u_cell_k(:));
        %         v_cell_k = v_cell_k - nanmedian(v_cell_k(:));
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
    
    if (plot_radial==1)
        save(cellvel_savename, 'x_cell', 'y_cell', 'u_cell','v_cell', 'ut', 'ur', "xc", "yc");
    else
        save(cellvel_savename, 'x_cell', 'y_cell', 'u_cell','v_cell', "xc", "yc");
    end
       % outputs: 
end