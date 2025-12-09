% Example code to use GetAllNearestPointsProfiles function

% add paths to MOPS and MOPS/toolbox to use the function
% eg. for a mac
% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% set a mop number
MopNumber=654%910;

%% load the 1m spatial avg "SA" survey data file 
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])

% print a list of the survey dates and sources to the command window
ListSA

% Use the nearest data point method to create profiles along the Mop
% transect

% set nearest point profile tolerances 
Ytol=50; % max dist (m) a survey point can be from the transect
Xtol=5; % max alongtransect gap (m) that will filled with linear interpolation

% GetNearestPointsProfile returns:
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1D(size(SA,2),numel(X1D)) = 2d matrix of nearest profiles for all SA surveys
%                                 no data = NaNs

[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

% make an example plots
idx=1:size(SA,2); % plot all surveys

col=jet(numel(idx));
m=0;
figure('position',[35          76        1329         714])
for n=idx
 m=m+1;
 pl(m)=plot(X1D,Z1D(n,:),'-','color',col(m,:),'linewidth',2,'DisplayName',...
        datestr(SA(n).Datenum,'mm/dd/yy'));hold on;
end
set(gca,'ylim',[-12 6],'fontsize',14);
xl=get(gca,'xlim');
plot(xl,xl*0+1.344,'k--');text(xl(2)-5,1.65,'MHW','fontsize',12);
set(gca,'xdir','reverse');
lg=legend(pl,'location','eastoutside','numcolumns',7,'fontsize',12);
title(lg,'Survey Dates')
box on;grid on;
title(['Mop ' num2str(MopNumber) ' Profiles']);
xlabel('Xshore Distance (m)');
ylabel('Elev (m, NAVD88)');

