% viewGMstruct.m
%load(GMf(GMmenu.Value).name);
load(fullfile(GMf(GMmenu.Value).folder,GMf(GMmenu.Value).name))

% T=struct2table(GM);
% sortedT = sortrows(T, 'Datenum');
% GM=table2struct(sortedT)';
% clear T sortedT;

% % Example to get indexes of Gps Jumbo surveys
% JumboIndexes=find(contains({GM.File}', 'jumbo','IgnoreCase',true)==1);
% if ~isempty(JumboIndexes)
% for n=JumboIndexes'
%     GM(n).Source='Gps-Jumbo';
% end
% end


figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[280 640 120 40],...
    'string',['MOP# ' num2str(GM.Mopnum) ' '],'foregroundcolor','b');

view2=uicontrol(fig2,'style','pushbutton','position',[400 650 70 30],...
    'string','1D','foregroundcolor','k','backgroundcolor','g',...
    'callback','delete(gca);delete(gca);delete(gca);;PlotCMstruct(GM,1);');
viewsd=uicontrol(fig2,'style','pushbutton','position',[480 650 70 30],...
    'string','2D','foregroundcolor','k','backgroundcolor','g',...
    'callback','delete(gca);delete(gca);delete(gca);PlotCMstruct(GM,2);');
view3=uicontrol(fig2,'style','pushbutton','position',[560 650 100 30],...
    'string','by Month','foregroundcolor','k','backgroundcolor','g',...
    'callback','delete(gca);delete(gca);delete(gca);PlotCMstruct(GM,2+DateMenu.Value);');

charmon = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}';

DateMenu=uicontrol(fig2,'style','popup','position',[670 640 90 40],...
   'string',charmon,'value',1,...
   'callback','cla;PlotCMstruct(GM,2+DateMenu.Value);');


cla;PlotCMstruct(GM,1);