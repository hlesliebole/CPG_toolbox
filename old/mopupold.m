%  Basic Viewer for Sand Survey Database NETCDF files

close all
toolboxpath=mfilename('fullpath');
toolfolder=fileparts(toolboxpath)
mopfolder=strrep(toolfolder,'toolbox','')

%global AreaTable AreaName Mop areamenu survmenu opendap ax1 ax2 H
%global fig fig2

% base path to the thredds beach survey catalog
% urlbase='http://thredds.cdip.ucsd.edu/thredds/catalog/test/sand/';
% % base opendap path to netcdf files
% opendap='http://thredds.cdip.ucsd.edu/thredds/dodsC/test/sand/';
% % path to MOP transect definition text file 
% mopurl='http://cdip.ucsd.edu/MOP_v1.1/CA_v1.1_transect_definitions.txt';

fig=figure('Name','CPGMOP Database','NumberTitle','off',...
    'position',[90 800 1300 60],'ToolBar','none','MenuBar','none');
set(fig,'Color',[223 206 157]/256);

fig2=figure('NumberTitle','off',...
     'position',[90 70 1300 680],'ToolBar','figure','MenuBar','none');

H=gobjects(4,1); % preallocate handle for profile map overlay pobjects;

%load_thredds_survey_netcdf_catalogs(urlbase)

%load_MOP_definitions(mopurl)

%launch_viewer_control_window

%------------------------------------------------------------------------
%function load_thredds_survey_netcdf_catalogs(urlbase)
%------------------------------------------------------------------------

% %global AreaTable AreaName 
% 
% % url of thredds sand survey files main catalog of survey area
% %  subfolders
% urlcat=[urlbase 'catalog.html'];
% % download catalog html info
% html=urlread(urlcat);
% % remove html characters
% html=regexprep(html, '<[^>]*>', '');
% % parse out the names of the survey area subfolders in the thredds catalog
% areas=regexp(html,';\S*/','match');
% %areas=strip(areas,';');
% areas=strrep(areas,';','');
% %areas=strip(areas,'/');
% areas=strrep(areas,'/','');
% num_areas=size(areas,2);
% 
% %---------------------------------------------------------------------
% % For each survey area, make a matlab table of survey
% %  dates and survey types by parsing survey file names in
% %  thredds folders
% %---------------------------------------------------------------------
% 
% % create cell array of table names for the survey areas
% AreaName=cell(num_areas,1);
% AreaTable=cell(num_areas,1);
% 
% % loop through areas and make tables
% 
% for a = 1:num_areas
%     
%     AreaName{a}=areas{a}; % table variable name
%     urldir=[urlbase areas{a} '/catalog.html']; % thredds file list url
%     html=urlread(urldir); % read in file list catalog as html
%     % remove html characters
%     html=regexprep(html, '<[^>]*>', '');
%     % parse out netcdf file names starting with sio
%     files=regexp(html,'sio\S*\.nc','match');
%     num_files=size(files,2);
%     % create cell arrays to store survey info for table
%     SurveyDate=cell(num_files,1);
%     SurveyType=cell(num_files,1);
%     SurveyGroup=cell(num_files,1);
%     SurveyName=cell(num_files,1);
%     % loop through files and add survey date, type and group to table cell
%     %   arrays
%     for f=1:num_files
%         info=regexprep(files{f},'\.nc', ''); %parse survey info from filename
%         info=regexprep(info,'_', ' ');
%         info=regexp(info,'\S*','match');
%         SurveyDate{f}=info{2};
%         SurveyType{f}=info{4};
%         SurveyGroup{f}=info{1};
%         SurveyName{f}=info{3};
%     end
%     
%     % create table for the survey area
%     AreaTable{a}=table(SurveyDate,SurveyType,SurveyName,SurveyGroup);
%     
% end
% 
% %end
% 
% %------------------------------------------------------------------------
% %function load_MOP_definitions(mopurl)
% %------------------------------------------------------------------------
% 
% %global Mop 
% 
% mopdef=webread(mopurl);
% i=strfind(mopdef,'--');   % find beginning of mop transect table
% mopdef=regexp(mopdef(i(end)+2:end),'\S*','match');
% mopdef=reshape(mopdef,[9 11594]);
% Name=mopdef(2,:)';
% BackLon=str2double(mopdef(3,:))';
% BackLat=str2double(mopdef(4,:))';
% OffLon=str2double(mopdef(5,:))';
% OffLat=str2double(mopdef(6,:))';
% Depth=str2double(mopdef(7,:))';
% Normal=str2double(mopdef(8,:))';
% Complex=str2double(mopdef(9,:))';
% Mop=table(Name,BackLon,BackLat,OffLon,OffLat,...
%    Depth,Normal,Complex); 

%end

%------------------------------------------------------------------------
%function launch_viewer_control_window
%------------------------------------------------------------------------

%global AreaTable AreaName areamenu survmenu mopmenu plotit 
%global aprev sprev fig 


AreaNum=1;
aprev=0;
sprev=0;

set(0,'DefaultUicontrolFontSize',20);
set(0,'DefaultUicontrolFontWeight','normal');

bc=([223 206 157]/256);
SAhead=uicontrol(fig,'style','text','position',[80 20 80 40],...
    'string','1m AVG ','foregroundcolor','b','backgroundcolor',bc);
SAf=dir([mopfolder 'M*SA.mat']);
SAmenu=uicontrol(fig,'style','popup','position',[10 1 190 40],...
   'string',vertcat(SAf.name));
viewSA=uicontrol(fig,'style','pushbutton','position',[200 10 50 30],...
    'string','view','foregroundcolor','k','backgroundcolor','g',...
    'callback','viewSAstruct');

% 
SGhead=uicontrol(fig,'style','text','position',[310 20 120 40],...
    'string','1m GRID','foregroundcolor','b','backgroundcolor',bc);
SGf=dir([mopfolder 'M*SG.mat']);
SGmenu=uicontrol(fig,'style','popup','position',[270 1 190 40],...
   'string',vertcat(SGf.name));
plotSG=uicontrol(fig,'style','pushbutton','position',[460 10 50 30],...
    'string','view','foregroundcolor','k','backgroundcolor','g',...
    'callback','viewSGstruct');

SMhead=uicontrol(fig,'style','text','position',[550 20 150 40],...
    'string','Morpho Params',...
    'foregroundcolor','b','backgroundcolor',bc);
SMf=dir([mopfolder 'M*SM.mat']);
SMmenu=uicontrol(fig,'style','popup','position',[530 1 195 40],...
   'string',vertcat(SMf.name));
plotSM=uicontrol(fig,'style','pushbutton','position',[720 10 50 30],...
    'string','view','foregroundcolor','k','backgroundcolor','g',...
    'callback','viewSMstruct');

GMhead=uicontrol(fig,'style','text','position',[790 20 150 40],...
    'string','Global Morpho ',...
    'foregroundcolor','b','backgroundcolor',bc);
GMf=dir([mopfolder 'M*GM.mat']);
GMmenu=uicontrol(fig,'style','popup','position',[780 1 195 40],...
   'string',vertcat(GMf.name));
plotGM=uicontrol(fig,'style','pushbutton','position',[970 10 50 30],...
    'string','view','foregroundcolor','k','backgroundcolor','g',...
    'callback','viewGMstruct');

SWhead=uicontrol(fig,'style','text','position',[1040 20 150 40],...
    'string','Waves',...
    'foregroundcolor','b','backgroundcolor',bc);
SWf=dir([mopfolder 'M*SW.mat']);
SWmenu=uicontrol(fig,'style','popup','position',[1030 1 195 40],...
   'string',vertcat(SWf.name));
plotSW=uicontrol(fig,'style','pushbutton','position',[1230 10 50 30],...
    'string','view','foregroundcolor','k','backgroundcolor','g',...
    'callback','viewSWstruct');

% 
% survmenu=uicontrol(fig,'style','popup','position',[250 1 500 40],...
%    'string',strcat(AreaTable{AreaNum}{:,1},{' - '},...
%     AreaTable{AreaNum}{:,2},{' - '},AreaTable{AreaNum}{:,3},{' - '},... 
%     AreaTable{AreaNum}{:,4}));
% 
% areamenu=uicontrol(fig,'style','popup','position',[10 1 230 40],...
%    'string',upper(AreaName),'callback',...
%    ['survmenu.Value=1;'...
%     'survmenu.String=join(AreaTable{areamenu.Value}{:,:}," : ");'...
%     'mopmenu.Value=1;mopmenu.String=''TBD''']);

% plotit=uicontrol(fig,'style','pushbutton','position',[1100 5 150 50],...
%     'string','Plot It','foregroundcolor','k','backgroundcolor','g',...
%     'callback','plot_survey');

%end





