%function [SRVhead,SRVmenu]=updateMOPsrv(fig,bc,SAf,MOPmenu)

% load first MOP SA file to get initial set of survey dates for menu
load(fullfile(SAf(MOPmenu.Value).folder,SAf(MOPmenu.Value).name));

% get indexes of Gps Jumbo surveys and add to string label
JumboIndexes=find(contains({SA.File}', 'jumbo','IgnoreCase',true)==1);
if ~isempty(JumboIndexes)
for n=JumboIndexes'
    SA(n).Source='Gps-Jumbo';
end
end

% SRVhead=uicontrol(fig2,'style','text','position',[250 20 180 40],...
%     'string','Survey Date','foregroundcolor','b','backgroundcolor',bc);
SRVmenu=uicontrol(fig2,'style','popup','position',[185 640 300 40],...
     'string',strcat(datestr(vertcat(SA.Datenum)),...
   repmat({' '},size(SA,2),1),char(SA.Source)),...
   'callback','viewStruct');

%end