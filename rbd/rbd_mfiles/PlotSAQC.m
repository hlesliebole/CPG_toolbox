%
% Makes a color dot plot of SA struct array elevationon data

%ScatterPlotBeachUTM(MobileLidarXunq,MobileLidarYunq,Zarea,'3d');
if exist('nav','var');delete(nav);end
if exist('ax1','var');delete(ax1);end

% ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')));
% fprintf(' 1. %i total Mobile LiDAR surveys\n',numel(ldx));
% fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
% fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));
sn=CurrentSurveyNumber;
x=vertcat(SA(ldx(sn)).X);
y=vertcat(SA(ldx(sn)).Y);
z=vertcat(SA(ldx(sn)).Z);

ax1=axes('position',[.1 .1 .8 .8]);
% plot colored points
minz=min(z);%quantile(z,.05);
maxz=max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
if dview == '3'
scp=scatter3(x,y,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 10 1]);
zoom off
elseif dview == 'T'
scp=scatter(x,y,12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 10 1]);
zoom on
elseif dview == 'N'
scp=scatter(x,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 1 1]);
zoom on
elseif dview == 'E'
scp=scatter(y,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 1 1],'xdir','reverse');
zoom on
end

cb=colorbar;cb.Label.String='Elevation';
set(ax1,'clim',[minz maxz],'fontsize',16);colormap(cm);
title(['MOP ' num2str(CurrentMopNumber) ' ' datestr(SA(ldx(sn)).Datenum) ' ' SA(ldx(sn)).Source] )
set(ax1,'color',[.5 .5 .5])

hold on;
    if exist('dp','var');delete(dp);end
if dview == '3'
    dp=plot3([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Y],...
     [QC(ldx(CurrentSurveyNumber)).Z],'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'T'
    dp=plot([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Y],...
     'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'N'
    dp=plot([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Z],...
     'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'E'
    dp=plot([QC(ldx(CurrentSurveyNumber)).Y],[QC(ldx(CurrentSurveyNumber)).Z],...
     'mx','MarkerSize',10,'linewidth',2);
end

set(ax1,'color',[.5 .5 .5])

% add navigation buttons to move bettwen surveys at a mop
SurveyNavigation;
