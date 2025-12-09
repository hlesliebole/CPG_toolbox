% Example code to use GetAllNearestPointsProfiles function

% add paths to MOPS and MOPS/toolbox to use the function

% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% load the desired Mop 1m spatial average SA struct array
load('M00511SA.mat','SA')

% set nearest point profile tolerances
%  For Lidar data, these can usually be set very tight
Ytol=3; % max dist (m) a survey point can be from the transect
Xtol=3; % max alongtransect gap (m) that will filled with linear interpolation

% GetAllNearestPointsProfiles returns:
%
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1D(size(SA,2),numel(X1D)) = 2d matrix of nearest profiles for all SA surveys
%                                 no data = NaNs

[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

% isolate the trk lidar data
Trkidx=find(contains({SA.Source}','Trk'));

% make a 3d figure of trk profiles
figure;hold on
for n=Trkidx'
plot3(X1D,SA(n).Datenum+X1D*0,Z1D(n,:));
end
plot3(X1D,SA(n).Datenum+X1D*0,Z1D(n,:),'m-','linewidth',2);
set(gca,'xdir','reverse','view',[ -29.1000   17.2105]);grid on
datetick('y')
