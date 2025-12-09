
% if in 3d view, switch to top mode to find points with ginput
if dview == '3'
  dview = 'T';
  PlotSAQC
end

% remove survey nav buttons
if exist('nav','var');delete(nav);end

% mfile for selecting a bathymetry point to delete

zoom off

% make dummy button for quitting

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','delete(hpop);delete(htxt);PlotSAQC;edit_menuSA','String','Done'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 670 800 25],'String', ...
'Click on individual points to plot their Z history. Click on Done button when finished.'); 
set(htxt,'ForegroundColor','red','fontsize',18);

xy=ginput(1);

% v=get(gca,'Ylim');
% if xy(2) < v(1)       
% zoom on
% return
% end

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

TimePlotFig=figure('position',[100 100 600 400]);
ax2=axes(TimePlotFig,'position',[.1 .1 .8 .8]);
% plot time series of point x(I), y(I)

% loop around

delete(hpop)
delete(htxt)
DeletePointsSA
