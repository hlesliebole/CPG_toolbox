% Example code to calcualate a MHW beach width for all the surveys at a
%  a single Mop transect using the CPG MOP M*SA.mat files

% add paths to MOPS and MOPS/toolbox to use the function
% eg. for a mac
% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts
clearvars
% set mop number, eg. cortez st in IB
Nmop=45;
Nmop=625;

% set nearest point profile tolerances 
Ytol=50; % max dist (m) a survey point can be from the transect
Xtol=5; % max alongtransect gap (m) that will filled with linear interpolation

% load the SA file of all the surveys
load(['M' num2str(Nmop,'%5.5i') 'SA.mat'],'SA');

% get nearest point profiles for all the sur ey sin the SA file
[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

%% loop through profiles and find MHW intersection with X1D as the
% definition of beach width

nsurv=size(Z1D,1);

for n=1:nsurv
 igood=find(~isnan(Z1D(n,:))); % profile points with valida data

 if numel(igood) > 1 % need at least 2 good profile points

 % find x intersections with MHW
 xMHW=intersections([X1D(1) X1D(end)],[1.344 1.344],X1D(igood),Z1D(n,igood));
 
   if ~isempty(xMHW)
     % keep the nearest and furthest intersection from the back beach
     % (often the same but useful for catching weird profiles
     BeachWidthMin(n)=xMHW(1); 
     BeachWidthMax(n)=xMHW(end);
   else
     BeachWidthMin(n)=NaN;
     BeachWidthMax(n)=NaN;
   end

 else
     BeachWidthMin(n)=NaN;
     BeachWidthMax(n)=NaN;
 end

end

% make a time series plot of the beach widths
Sdatetime=datetime([SA.Datenum],'convertfrom','datenum');

figure('position',[ 285   195   858   414]); hold on;
plot(Sdatetime,BeachWidthMin,'+-','DisplayName','Min MHW Beach Width');
plot(Sdatetime,BeachWidthMax,'o--','DisplayName','Max MHW Beach Width');
legend('location','northoutside')
grid on;
xlabel('Date');ylabel('Distance from Mop Back Beach Point (m)')