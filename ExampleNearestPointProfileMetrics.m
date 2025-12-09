% Example code to use GetAllNearestPointsProfiles function
%   to generate MOP profiles from a CPG Mop SA.mat file and
%   calculate MHW beach width and volume metrics

% add paths to MOPS and MOPS/toolbox to use the function

% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% set Mop number
MopNumber=42;

% load the desired Mop 1m spatial average SA struct array
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA')

% set nearest point profile tolerances 
Ytol=50; % max dist (m) a survey point can be from the transect
Xtol=5; % max alongtransect gap (m) that will filled with linear interpolation

% GetNearestPointsProfile returns:
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1D(size(SA,2),numel(X1D)) = 2d matrix of nearest profiles for all SA surveys
%                                 no data = NaNs

[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

% Isolate profiles bettwen Nov 21 and May 24

idx=find([SA.Datenum] >= datenum(2021,11,1) & [SA.Datenum] < datenum(2024,5,1));

% make a 3d figure of jumbo profiles

figure('position',[225          60        1041         735]);hold on

subplot(3,1,1);hold on;
col=jet(numel(idx));
m=0;
for n=idx
 m=m+1;
%  MHW beach width from MOP Back Beach Point
xMHW(m)=intersections([X1D(1) X1D(end)],[1.344 1.344],X1D,Z1D(n,:));
% Profile sand volume above MHW
vdx=find(Z1D(n,:) > 1.344);
VolMHW(m)=sum(Z1D(n,vdx)-1.344,'omitnan');
% plot profile
pl(m)=plot(X1D,Z1D(n,:),'-','color',col(m,:),'linewidth',2,'DisplayName',...
        datestr(SA(n).Datenum,'mm/dd/yy'));hold on;
end
xl=get(gca,'xlim');
plot(xl,xl*0+1.344,'k--');text(xl(2)-5,1.55,'MHW')
set(gca,'xdir','reverse');
legend(pl,'location','eastoutside','numcolumns',3);
box on;grid on;
title(['Mop ' num2str(MopNumber) ' Profiles'])
xlabel('Xshore Distance (m)');
ylabel('Elev (m, NAVD88)')

% plot beach width time series
subplot(3,1,2);
plot(datetime([SA(idx).Datenum],'convertfrom','datenum'),xMHW,'*-');
box on;grid on;
xlabel('Survey Date');
ylabel('MHW Xshore Location, m)')
title('MHW Beach Width')

% plot beach volume time series
subplot(3,1,3);
plot(datetime([SA(idx).Datenum],'convertfrom','datenum'),VolMHW,'*-');
box on;grid on;
xlabel('Survey Date');
ylabel('Volume (m^{3}/m-shoreline')
title('Beach Volume Above MHW ')


