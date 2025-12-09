%
%  Main MATLAB script for viewing and editing .rdb files
%

% put subdirectory with the other mfiles in path

clearvars
close all

path(path,'./rbd_mfiles')

% get any default settings from defaults.m

rbd_defaults
close

% initialize stuff
ns=0;
nb=0;
nn=0;

dview='3';

format compact

% load first file

% try and use second monitor
% p = get(0, "MonitorPositions");
% if size(p,1) > 1
%     f.Position = p(2, :); % second display
%     f.WindowState = "maximized";
% end
% set(0, "DefaultFigurePosition", f.Position)

MainFig=figure('position',[250  103  1183 694]);

% default start is 654 fletcher cove
CurrentMopNumber=654;
CurrentSurveyNumber=1;

% get the SA and QC struct arrays for the starting Mop
GetMopMobileLidarSA


% SAmatfile=['M' num2str(CurrentMopNumber,'%5.5i') 'SA.mat'];
% fprintf('Loading SA struct array from %s\n',SAmatfile)
% load(SAmatfile);
% QCmatfile=['M' num2str(CurrentMopNumber,'%5.5i') 'QC.mat'];
% % get QC struct array if exist, otherwise make one to mirror SA
% %  but without any matching x,y,z qc removal points
% fprintf('Loading QC struct array from %s\n',QCmatfile)
% if exist("QCmatfile","file")
%     load(QCmatfile);
% else
%     fprintf('No QC file found. Making empty QC struct array.\n')
%     QC=SA;
%     for n=1:size(QC,2)
%         QC(n).X=[]; QC(n).Y=[]; QC(n).Z=[];
%     end
% end

% PlotSAQC
% edit_menuSA

% p = get(0, "MonitorPositions");
% f.Position = p(1, :); % first display
% f.WindowState = "maximized";
% set(0, "DefaultFigurePosition", "remove")


%close(gcf)
