% viewSGstruct.m
%load(SGf(SGmenu.Value).name);
load(fullfile(SGf(SGmenu.Value).folder,SGf(SGmenu.Value).name))

% resort struct array by date
T=struct2table(SG);
sortedT = sortrows(T, 'Datenum');
SG=table2struct(sortedT)';
clear T sortedT;

% Example to get indexes of Gps Jumbo surveys
JumboIndexes=find(contains({SG.File}', 'jumbo','IgnoreCase',true)==1);
if ~isempty(JumboIndexes)
for n=JumboIndexes'
    SG(n).Source='Gps-Jumbo';
end
end


figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
    'string',['MOP# ' num2str(SG(1).Mopnum) ' '],'foregroundcolor','b');
DateMenu=uicontrol(fig2,'style','popup','position',[120 640 280 40],...
   'string',strcat(datestr(vertcat(SG.Datenum)),...
   repmat({' '},size(SG,2),1),char(SG.Source)),...
   'callback','cla;PlotCMstruct(SG,DateMenu.Value);');
view2=uicontrol(fig2,'style','pushbutton','position',[400 650 50 30],...
    'string','Top','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(2)');
viewsd=uicontrol(fig2,'style','pushbutton','position',[460 650 70 30],...
    'string','Side','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(0,0)');
view3=uicontrol(fig2,'style','pushbutton','position',[540 650 50 30],...
    'string','3D','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(3)');

cla;PlotCMstruct(SG,DateMenu.Value);