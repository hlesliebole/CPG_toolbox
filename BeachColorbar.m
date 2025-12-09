
ax=gca;
% axpos=get(gca,'position'); % get current axes

% Get Beach colormap variables
load BeachColorMap.mat
% ic=find(round(10000*BeachColorMap(:,1)) == 5307);
% BeachColorMap(ic,:)=BeachColorMap(ic,:)*.7;
% ic=find(round(10000*BeachColorMap(:,1)) == 4507);
% BeachColorMap(ic,:)=BeachColorMap(ic,:)*.85;

% make Beach colormap colorbar
colormap(BeachColorMap);
BeachColorBar=colorbar('AxisLocation','in'); 
set(gca,'clim',BeachColorBarLims); 
BeachColorBar.Title.String=BeachColorBarTitleString;
BeachColorBar.Ticks=BeachColorBarTicks;
BeachColorBar.TickLabels=BeachColorBarTickLabels;
BeachColorBar.TickLength=0;

axpos=get(ax,'position'); % get plot axes position
set(ax,'position',axpos) % fix plot axes position

% shift colorbar right  to accomodate left hand labels
newpos=axpos(1)+axpos(3);newpos=newpos+0.45*(1-newpos);
pos=BeachColorBar.Position;pos(1)=newpos;BeachColorBar.Position=pos;
% 
% % make left y-axis with navd88 elevation labels
NavdAxes=axes('position',BeachColorBar.Position,'color','none','xtick',[],...
    'ytick',[-8:5],'ylim',[-8 5],'TickLength',[0 .01],'YAxisLocation',...
    'right');
% 
axes(ax); % return to plotting axes