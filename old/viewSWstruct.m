% viewSWstruct.m
%load(SWf(SWmenu.Value).name);
load(fullfile(SWf(SWmenu.Value).folder,SWf(SWmenu.Value).name))

T=struct2table(SW);
sortedT = sortrows(T, 'Datenum');
SW=table2struct(sortedT)';
clear T sortedT;

% % Example to get indexes of Gps Jumbo surveys
% JumboIndexes=find(contains({SW.File}', 'jumbo','IgnoreCase',true)==1);
% if ~isempty(JumboIndexes)
% for n=JumboIndexes'
%     SW(n).Source='Gps-Jumbo';
% end
% end


figure(fig2);clf;
% DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
%     'string',['MOP# ' num2str(SW(1).Mopnum) ' '],'foregroundcolor','b');
% DateMenu=uicontrol(fig2,'style','popup','position',[120 640 280 40],...
%    'string',strcat(datestr(vertcat(SW.Datenum)),...
%    repmat({' '},size(SW,2),1),char(SW.Source)),...
%    'callback','cla;PlotCMstruct(SW,DateMenu.Value);');
% view2=uicontrol(fig2,'style','pushbutton','position',[400 650 70 30],...
%     'string','1D','foregroundcolor','k','backgroundcolor','g',...
%     'callback','delete(gca);delete(gca);delete(gca);;PlotCMstruct(SW,1);');
% viewsd=uicontrol(fig2,'style','pushbutton','position',[480 650 70 30],...
%     'string','2D','foregroundcolor','k','backgroundcolor','g',...
%     'callback','delete(gca);delete(gca);delete(gca);PlotCMstruct(SW,2);');
% view3=uicontrol(fig2,'style','pushbutton','position',[540 650 50 30],...
%     'string','3D','foregroundcolor','k','backgroundcolor','g',...
%     'callback','view(3)');

cla;PlotCMstruct(SW);