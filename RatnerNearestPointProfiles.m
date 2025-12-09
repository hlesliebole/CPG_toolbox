% Example code to use GetAllNearestPointsProfiles function

% add paths to MOPS and MOPS/toolbox to use the function

% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% load the desired Mop 1m spatial average SA struct array
load(['M02672SA.mat'],'SA')

% set nearest point profile tolerances 
Ytol=30; % max dist (m) a survey point can be from the transect
Xtol=5; % max alongtransect gap (m) that will filled with linear interpolation

% GetNearestPointsProfile returns:
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1D(size(SA,2),numel(X1D)) = 2d matrix of nearest profiles for all SA surveys
%                                 no data = NaNs

[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

% isolate the jumbos
%Jidx=find(contains({SA.File}','jumbo','IgnoreCase',true)==1);
% everything
Jidx=(1:size(SA,2))';%find(contains({SA.File}','jumbo','IgnoreCase',true)==1);

% make a 3d figure of jumbo profiles
figure;hold on
for n=Jidx'
plot(X1D,Z1D(n,:));
end
set(gca,'xdir','reverse')%,'view',[ -29.1000   17.2105]);grid on
%datetick('y')
