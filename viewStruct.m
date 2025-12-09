
%-------------------------------------------------------------

if Smenu.Value == 1
    
load(fullfile(SAf(MOPmenu.Value).folder,SAf(MOPmenu.Value).name))

% T=struct2table(SA);
% sortedT = sortrows(T, 'Datenum');
% SA=table2struct(sortedT)';
% clear T sortedT;


figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
    'string',['MOP# ' num2str(SA(1).Mopnum) ' '],'foregroundcolor','b');
view2=uicontrol(fig2,'style','pushbutton','position',[540 650 50 30],...
    'string','Top','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(2)');
viewsd=uicontrol(fig2,'style','pushbutton','position',[610 650 70 30],...
    'string','Side','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(0,0)');
view3=uicontrol(fig2,'style','pushbutton','position',[690 650 50 30],...
    'string','3D','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(3)');

% get indexes of Gps Jumbo surveys and add to string label
JumboIndexes=find(contains({SA.File}', 'jumbo','IgnoreCase',true)==1);
if ~isempty(JumboIndexes)
for n=JumboIndexes'
    SA(n).Source='Gps-Jumbo';
end
end

SRVmenu=uicontrol(fig2,'style','popup','position',[185 640 300 40],...
     'string',strcat(datestr(vertcat(SA.Datenum)),...
   repmat({' '},size(SA,2),1),char(SA.Source)),'Value',Snum,...
   'callback','cla;Snum=SRVmenu.Value;PlotCMstruct(SA,Snum);');

Snum=SRVmenu.Value; % generic survey index number passed between plot figs
cla;PlotCMstruct(SA,Snum);

end

%-------------------------------------------------------------

if Smenu.Value == 2
    
load(fullfile(SGf(MOPmenu.Value).folder,SGf(MOPmenu.Value).name))

% T=struct2table(SG);
% sortedT = sortrows(T, 'Datenum');
% SG=table2struct(sortedT)';
% clear T sortedT;

figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
    'string',['MOP# ' num2str(SG(1).Mopnum) ' '],'foregroundcolor','b');
view2=uicontrol(fig2,'style','pushbutton','position',[540 650 50 30],...
    'string','Top','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(2)');
viewsd=uicontrol(fig2,'style','pushbutton','position',[610 650 70 30],...
    'string','Side','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(0,0)');
view3=uicontrol(fig2,'style','pushbutton','position',[690 650 50 30],...
    'string','3D','foregroundcolor','k','backgroundcolor','g',...
    'callback','view(3)');

GSRVmenu=uicontrol(fig2,'style','popup','position',[185 640 300 40],...
     'string',strcat(datestr(vertcat(SA.Datenum)),...
   repmat({' '},size(SA,2),1),char(SA.Source)),...
   'Value',Snum,...
   'callback','cla;Snum=GSRVmenu.Value;PlotCMstruct(SG,Snum);');

   
cla;PlotCMstruct(SG,Snum);

end

%-------------------------------------------------------------
if Smenu.Value == 3
    
load(fullfile(SMf(MOPmenu.Value).folder,SMf(MOPmenu.Value).name))

% T=struct2table(SM);
% sortedT = sortrows(T, 'Datenum');
% SM=table2struct(sortedT)';
% clear T sortedT;

figure(fig2);clf;
DateHead=uicontrol(fig2,'style','text','position',[10 640 120 40],...
    'string',['MOP# ' num2str(SM(1).Mopnum) ' '],'foregroundcolor','b');
view2=uicontrol(fig2,'style','pushbutton','position',[120 650 70 30],...
    'string','Prev','foregroundcolor','k','backgroundcolor','g',...
    'callback','if MSRVmenu.Value > 1;MSRVmenu.Value=MSRVmenu.Value-1;cla;Snum=MSRVmenu.Value;PlotCMstruct(SM,Snum);end');
viewsd=uicontrol(fig2,'style','pushbutton','position',[480 650 70 30],...
    'string','Next','foregroundcolor','k','backgroundcolor','g',...
    'callback','if MSRVmenu.Value < size(SM,2);MSRVmenu.Value=MSRVmenu.Value+1;cla;Snum=MSRVmenu.Value;PlotCMstruct(SM,Snum);end');

MSRVmenu=uicontrol(fig2,'style','popup','position',[185 640 300 40],...
     'string',strcat(datestr(vertcat(SA.Datenum)),...
   repmat({' '},size(SA,2),1),char(SA.Source)),'Value',Snum,...
   'callback','cla;Snum=MSRVmenu.Value;PlotCMstruct(SM,Snum);');

cla;PlotCMstruct(SM,Snum);

end

%-------------------------------------------------------------
if Smenu.Value == 4
    
load(fullfile(GMf(MOPmenu.Value).folder,GMf(MOPmenu.Value).name))

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

end


%-------------------------------------------------------------
if Smenu.Value == 5
    
load(fullfile(SWf(MOPmenu.Value).folder,SWf(MOPmenu.Value).name))

T=struct2table(SW);
sortedT = sortrows(T, 'Datenum');
SW=table2struct(sortedT)';
clear T sortedT;

figure(fig2);clf;

cla;PlotCMstruct(SW);

end