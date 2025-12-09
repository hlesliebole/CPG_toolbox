% MobileLidarXutm=vertcat(SA(ldx).X);
% MobileLidarYutm=vertcat(SA(ldx).Y);
% MobileLidarZnavd88=vertcat(SA(ldx).Z);
% % make a survey date to go with every survey point
% for n=1:size(SA,2)
%     SA(n).Dates=repmat(SA(n).Datenum,size(SA(n).X));
%     %SA(n).Dates=repmat(n,size(SA(n).X));
% end
% MobileLidarDatenum=vertcat(SA(ldx).Dates);
% %figure;histogram(MobileLidarZnavd88)
% 
% [uxy,ia,ic]=unique([MobileLidarXutm,MobileLidarYutm],'rows');
% fprintf(' 4. %i unique LiDAR Xutm,Yutm points in area.\n',numel(ia));
% 
% %% ------------------------------------------------------------------
% 
% 
% % remove survey nav buttons
% if exist('nav','var');delete(nav);end
% 
% % mfile for selecting a bathymetry point to delete
% 
% zoom off
% 
% % make dummy button for quitting
% 
% hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
% 'Callback','delete(hpop);delete(htxt);PlotSAQC;edit_menuSA','String','Done'); 
% 
% htxt=uicontrol(gcf,'Style','text','Position',...
% [150 670 800 25],'String', ...
% 'Click on individual points to plot their entire Z time series. Click on Done button when finished.'); 
% set(htxt,'ForegroundColor','red','fontsize',18);
% 

%% Find all the data points at the mouse input point


% find closet lat lon point
%vw=get(ax1,'view');
% looking from top, match with x and y
if dview == 'T'
 d=sqrt((x-xy(1)).^2+(y-xy(2)).^2);
 [Y,I]=min(d);
end
% looking north, match with x and z
if dview == 'N'
 d=sqrt((x-xy(1)).^2+(z-xy(2)).^2);
 [Y,I]=min(d);
end
% looking east, match with y and z
if dview == 'E'
 d=sqrt((y-xy(1)).^2+(z-xy(2)).^2);
 [Y,I]=min(d);
end

% find a plot a z data for this point

ipts=find(ismember([MobileLidarXutm';MobileLidarYutm']',[x(I)';y(I)']','rows'));

TimeFig=figure('position',[250 107 1178 371]);
ax2=axes('position',[.15 .1 .8 .8]);set(gca,'fontsize',14);ylabel('elevation (m,navd88)')
plot(datetime(MobileLidarDatenum(ipts),'convertfrom','datenum'),MobileLidarZnavd88(ipts),'*-');
hold on;
isrv=find(MobileLidarDatenum(ipts) == SA(ldx(CurrentSurveyNumber)).Datenum);
plot(datetime(MobileLidarDatenum(ipts(isrv)),'convertfrom','datenum'),...
    MobileLidarZnavd88(ipts(isrv)),'mo','markersize',10,'linewidth',2);

grid on;

hpop2=uicontrol(TimeFig,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','delete(hpop2);delete(ax2);delete(TimeFig);figure(MainFig);SelectZtimeHistorySA','String','Clear'); 