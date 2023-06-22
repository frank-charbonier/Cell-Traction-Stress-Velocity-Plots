function plot_cellvel(cellname, domainname, cellvel_savename, pix_size, time_increment, max_vel, plot_radial)
arguments
    % Name of multipage tif file to plot cells.
    % Set to 'none' if there is no cell image
    cellname = 'cells.tif';
    % Name of domain. This is where cells are located. Set to [] if no domain
    domainname = 'domain.tif';
    % Name of cell velocity data to load
    cellvel_savename = 'cellvel_processed.mat';
    % Pixel size [microns]
    pix_size = 1.3;
    % Time between images
    time_increment = 10; % min
    % Max velocity for color plots
    max_vel = 0.25;   % units: um/min
    % Assay format
    % Set to 0 to plot x and y components, otherwise plot radial and tangential
    plot_radial = 0;
end
% PLOT_CELLVEL Plot cell velocities
%
% First run compute_cellvel.m, then use this for plotting
%
% Written by Jacob Notbohm, Univerity of Wisconsin-Madison, 2015-2020
% Adapted by Frank Charbonier, Stanford University, 2023

% clear;
close all;
% clc;

%% --- USER INPUTS ---
% Typically, the quiver plot, the number of quivers needs to be reduced so
% that they can be seen more clearly. Downsample the number of data points
% for plotting quivers by this factor
qd = 8;
% Quiver size
quiver_size=3;
% Header of name to save plots. This can contain a directory listing
dirname = 'cell_velocity'; % Name of a folder to put plots in
% dirname = 'cell_velocity_qd_8_quiver_size_3'; % Name of a folder to put plots in
savenameheader = [dirname,'/t_']; % Header of file name to save
% Set to [] to make figure visible. Set to 1 to make figure invisible.
invisible = 1;
% Name of file to save outputs
savename_data = 'cellvel_processed.mat';

%% --- MAKE PLOTS ---

% Make directory to save plots
if exist(dirname,'dir')==7
    rmdir(dirname,'s');
end
mkdir(dirname);

% Load data
load(cellvel_savename);

% Number of correlations
K = size(u_cell,3);

cmap = load('map_cold-hot.dat'); % Load hot-cold colormap data

for k=1:K
    hf = make_fig([0.2 0.2 2 1.4]);
    set(hf,'DefaultAxesPosition',[0.03 0.05 .95 .89]);
    if invisible
        set(hf,'visible','off');
    end
    
    % Cell images
    if strcmp(cellname,'none')==0
        % First cell image
        subplot(2,3,1)
        im_k = imread(cellname,k);
        [M, N] = size(im_k);
        imagesc([0 N]*pix_size,[0 M]*pix_size,im_k); colormap(gca,'gray');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        if pix_size == 1
            xlabel('pix'); ylabel('pix');
        else
            xlabel('\mum'); ylabel('\mum');
        end
        title('Reference')
        % Second cell image
        subplot(2,3,4)
        im_k = imread(cellname,k+1);
        [M, N] = size(im_k);
        imagesc([0 N]*pix_size,[0 M]*pix_size,im_k); colormap(gca,'gray');
        axis xy; axis equal; axis tight; set(gca,'box','off');
        xlabel('\mum'); ylabel('\mum');
        title('Current')
    end
    
    % Load cell velocity for current frame
    u_cell_k = u_cell(:,:,k);
    v_cell_k = v_cell(:,:,k);
    
    % Plot speed
    subplot(2,3,5)
    u_cell_mag = sqrt(u_cell_k.^2 + v_cell_k.^2);
    % imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_mag);
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_mag, ...
        "AlphaData",~isnan(u_cell_mag)); % use AlphaData property to make naan values transparent
    hold on
    x_cell2=downsample(x_cell,qd); x_cell2=downsample(x_cell2',qd)';
    y_cell2=downsample(y_cell,qd); y_cell2=downsample(y_cell2',qd)';
    u_cell_k2=downsample(u_cell_k,qd); u_cell_k2=downsample(u_cell_k2',qd)';
    v_cell_k2=downsample(v_cell_k,qd); v_cell_k2=downsample(v_cell_k2',qd)';
    quiver(x_cell2,y_cell2,u_cell_k2,v_cell_k2, quiver_size,... % The scalar number is the relative scaling (length) of the quivers
        'color', [0 0 0],'linewidth',1); %You may have to adjust the color of the quivers to show up agaist the colormap
    caxis([0 max_vel]);colormap(gca, brewermap([],'YlOrRd')); colorbar;
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('Cell speed');
    
    % Plot x velocity
    subplot(2,3,2)
    % imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_k);
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],u_cell_k, ...
        "AlphaData",~isnan(u_cell_k)); % use AlphaData property to make naan values transparent
    caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
    set(gca, 'Color', 'k')   % Set plot background to black
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('X velocity');
    
    % Plot y velocity
    subplot(2,3,3)
    % imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],v_cell_k);
    imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],v_cell_k, ...
        "AlphaData",~isnan(u_cell_k)); % use AlphaData property to make naan values transparent
    caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
    set(gca, 'Color', 'k')   % Set plot background to black
    axis xy; axis equal; axis tight; set(gca,'box','off');
    if pix_size == 1
        xlabel('pix'); ylabel('pix');
    else
        xlabel('\mum'); ylabel('\mum');
    end
    title('Y velocity');
    
    % OPTIONAL: radial and angular compoents of velocity instead of x and y
    if (plot_radial==1)
        % Load cell velocity for current frame
        ur_k = ur(:,:,k);
        ut_k = ut(:,:,k);
        % Radial velocity
        subplot(2,3,2);
    
        % r_grid = sqrt( (x-xc).^2 + (y-yc).^2 ); % Gridpoints
        % theta = atan2( (y_cell-yc), (x_cell-xc) );
    
%         imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ur_k);
%         caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
        imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ur_k, ...
            "AlphaData",~isnan(u_cell_k)); % use AlphaData property to make naan values transparent
        caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
        set(gca, 'Color', 'k')   % Set plot background to black
        axis xy; axis equal; axis tight; set(gca,'box','off');
        xlabel('\mum','fontsize',11); ylabel('\mum','fontsize',11);
        title('Cell radial velocity');
    
        % Angular velocity
        subplot(2,3,3)    
%         imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ut_k);
%         caxis([-max_vel max_vel]); colormap(gca,cmap); colorbar;
        imagesc([min(x_cell(:)) max(x_cell(:))],[min(y_cell(:)) max(y_cell(:))],ut_k, ...
            "AlphaData",~isnan(u_cell_k)); % use AlphaData property to make naan values transparent
        caxis([-max_vel max_vel]); colormap(gca, cmap); colorbar;
        set(gca, 'Color', 'k')   % Set plot background to black
        axis xy; axis equal; axis tight; set(gca,'box','off');
        xlabel('\mum','fontsize',11); ylabel('\mum','fontsize',11);
        title('Cell angular velocity');

    end
    % Save
     set(gcf,'PaperPositionMode','auto','InvertHardCopy','off');
    print('-dpng','-r300',[savenameheader,num2str(k,'%0.3d'),'-',num2str(k+1,'%0.3d')]);
    
    close(hf);
    
end