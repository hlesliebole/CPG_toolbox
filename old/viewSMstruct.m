% viewSMstruct.m
%load(SMf(SMmenu.Value).name);
load(fullfile(SMf(SMmenu.Value).folder,SMf(SMmenu.Value).name))

T=struct2table(SM);
sortedT = sortrows(T, 'Datenum');
SM=table2struct(sortedT)';
clear T sortedT;

% Example to get indexes of Gps Jumbo surveys
JumboIndexes=find(contains({SM.File}', 'jumbo','IgnoreCase',true)==1);
if ~isempty(JumboIndexes)
for n=JumboIndexes'
    SM(n).Source='Gps-Jumbo';
end
end


figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
    'string',['MOP# ' num2str(SM(1).Mopnum) ' '],'foregroundcolor','b');
DateMenu=uicontrol(fig2,'style','popup','position',[120 640 280 40],...
   'string',strcat(datestr(vertcat(SM.Datenum)),...
   repmat({' '},size(SM,2),1),char(SM.Source)),...
   'callback','cla;PlotCMstruct(SM,DateMenu.Value);');
view2=uicontrol(fig2,'style','pushbutton','position',[400 650 70 30],...
    'string','Prev','foregroundcolor','k','backgroundcolor','g',...
    'callback','if DateMenu.Value > 1;DateMenu.Value=DateMenu.Value-1;cla;PlotCMstruct(SM,DateMenu.Value);end');
viewsd=uicontrol(fig2,'style','pushbutton','position',[480 650 70 30],...
    'string','Next','foregroundcolor','k','backgroundcolor','g',...
    'callback','if DateMenu.Value < size(SM,2);DateMenu.Value=DateMenu.Value+1;cla;PlotCMstruct(SM,DateMenu.Value);end');
% view3=uicontrol(fig2,'style','pushbutton','position',[540 650 50 30],...
%     'string','3D','foregroundcolor','k','backgroundcolor','g',...
%     'callback','view(3)');

cla;PlotCMstruct(SM,DateMenu.Value);