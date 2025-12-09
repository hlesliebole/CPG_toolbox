% Example matlab code to load previously processed SHOALS bottom 
%   roughness and elevation data and plot on a google map figure.
%
%  uses the additional m-scripts:
%
%       utm2deg.m
%       plot_google_map.m
%
%   and .mat shoals data files 
%
%       M00660to00723RC.mat
%       M00474to00495RC.mat 
% 

%%------------------------------------------------------------
%  Load a processed SHOALS data file for the desired Mop range.
%
%   The .mat files contain a matlab struct array called "RC"
%   for "R"ougness "C"omposite data. The array has the following
%   elements:
%      RC.Mopnum
%      RC.SpatialWindowSize;
%      RC.SpatialMinPoints;
%      RC.Xutm 
%      RC.Yutm  
%      RC.Z % spatial filter avg elev NAVD88
%      RC.Sigma   % spatial filter standard deviation "roughness"
%
%    where RC has N entries, one for each Mop area in the struct array.
%     eg. RC(1).Sigma is the roughnees values for the shoals points in 
%          the first Mop area etc.
%%------------------------------------------------------------

% Option 1 for the Pt La Jolla area, from the Children's pool on
%  the west side of the point (Mop 660) to the Marine Room restaurant 
%  on the south end of LJ shores bech (Mop 723).
%
%load M00474to00495RC.mat 

% Option 2 for Carslbad SB and San Elijo SB area, from seaside reef at
%  the south end of Carlsbad SB (Mop 660) to Moonligh SB in Encinitas (Mop723).
%
addpath /Users/William/Desktop/Connor
load M00660to00723RC.mat

% grandview
%load M00755to00759RC.mat

% Mop data is normally saved as Xutm and Yutm (UTM eastings and northings)
%  coordinates in the UTM "11 S" zone, which is the most convenient for
%  calculating areas and volumes because UTM is a 1m x 1m rectangular coord
%  system. But here, convert to x,y (Longitude,Latitude)
%  coordinates to plot them on a google map.  Using "vertcat" with the RC
%  struct array element combines the values for ALL the Mops in the array.
[y,x]=utm2deg(vertcat(RC.Xutm),vertcat(RC.Yutm),...
    repmat('11 S',[length(vertcat(RC.Xutm)) 1]));

%-----------------------------------
% Plot roughness values on google map
%-----------------------------------
z=vertcat(RC.Sigma); % use roughness for z
figure('position',[110 503 713 855]);
ax1=axes;set(ax1,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);hold on;
% plot google map
plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)
% overlay roughness points as colored scatter points
zmin=0;zmax=0.5;set(gca,'clim',[zmin  zmax]); % show 0-0.5 stdev roughness range
zrange=zmax-zmin; % set max range of coloring
% do some scaling stuff to assign colors to roughness values
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % standard plotting check for any NaN no data points 
cm =jet(64);colormap(jet(64)); % use the matlab jet colormap
% overlay scatter plot points on map
scatter(x(idx), y(idx), 1, cm(ceil(zscaled(idx)),:), 'filled');
set(gca,'dataaspectratio',[1 cosd(y(1)) 1]);
set(ax1,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
cb=colorbar('horizontal');box on;
cb.Label.String='2D Spatial Standard Deviation "Substrate Roughness" (m)';
xlabel('Longitude');ylabel('Latitude');
title([{'USACE SHOALS Topo-Bathymetric LiDAR Elevation Standard Deviations'},...
    {'Using a 11m x 11m 2D Moving Filter'}]);

%-----------------------------------
% Plot mean depth values on a google map
%-----------------------------------
z=vertcat(RC.Z); % use mean filtered window depths as z
% only show below MSL
% convert to MSL
z=z-0.774;
x(z > 0)=[];y(z > 0)=[];z(z > 0)=[]; % toss data > msl
x(z < -10)=[];y(z < -10)=[];z(z < -10)=[]; % toss data < -10m
%x(z > 0.774)=[];y(z > 0.774)=[];z(z > 0.774)=[]; % toss data > msl
figure('position',[523    97   319   689]);
%ax1=axes;set(ax1,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);hold on;
ax1=axes;set(ax1,'xlim',[-117.2886 -117.2772],'ylim',[32.9970   33.0167]);hold on;

% plot google map
plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)
% overlay elevation points
zmin=-10;zmax=0.;set(gca,'clim',[zmin zmax]); % show -20 to MSL range
zrange=zmax-zmin; % set max range of coloring
% do some scaling stuff to assign colors to roughness values
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % standard plotting check for any NaN no data points 
cm =jet(64);colormap(jet(64)); % use the matlab jet colormap
% overlay scatter plot points on map
scatter(x(idx), y(idx), 2, cm(ceil(zscaled(idx)),:), 'filled');
set(gca,'dataaspectratio',[1 cosd(y(1)) 1]);
set(ax1,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
cb=colorbar('horizontal');box on;
cb.Label.String='Elevation (m, MSL)';
xlabel('Longitude');ylabel('Latitude');
% title([{'USACE SHOALS Topo-Bathymetric LiDAR Elevations'},...
%     {'Using a 11m x 11m 2D Moving Filter'}]);
set(ax1,'xlim',[-117.2886 -117.2772],'ylim',[32.9970   33.0167]);
