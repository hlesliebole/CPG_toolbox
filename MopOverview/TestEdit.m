
clearvars

%% find and add some paths to the CPG MOP data and toolbox directories 
% Use the fletcher cove M00654SA.mat file location to set paths
CpgMopFilePath = fileparts(which('M00654SA.mat'));
if ~exist('CpgMopFilePath','var')
    fprintf('Path to Fletcher M00654SA.mat not found. Stopping.')
    return
else
    % add to path
    addpath(CpgMopFilePath)
    % add MOPS/toolbox to path 
    addpath([CpgMopFilePath filesep 'toolbox']);
    addpath([CpgMopFilePath filesep 'toolbox' filesep 'MopOverview']);
    addpath([CpgMopFilePath filesep 'toolbox' filesep 'MopSAeditor']);
  
end


%% Default starting Mop number for volume editor
CurrentMopNumber=654;
VolumeElevCutoff=0.774; % MSL subaerial volume
MinFractionSurveysWithData=0.05;
MinFractionalCoverageAboveElevCutoff=0.5;

VolFig=figure('position',[5          72        1426         717]);

%% make some editable gui fields to define volume estimates
%
%   1) Mop number (CurrentMopNumber)
%
%   2) Elevation to use for defining subaerial volumes (VolumeElevCutoff)
%
%                   Sparse Data Tolerances
%
%   3) Minimum allowable fraction of the surveys where a Mop x,y point was 
%       observed with lidar.  This defines the fixed max control area for 
%       all the volume estimates. (MinFractionSurveysWithData)
%
%   4) For individual surveys, Minimum allowable fraction of 1m spatial 
%       resolution x,y points landward of the elevation in 2) that have 
%       valid lidar data. (MinFractionalCoverageAboveElevCutoff)
%       

%%  area sufficient data to make a volume estimate.

%%  Change Mop to Edit

Vnav(1)=uicontrol(VolFig,'style','text','position',[10 690 70 25],...
    'string','MOP #','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',18);

Vnav(2)=uicontrol(VolFig,'style','edit','position',[85 690 100 25],...
    'string',CurrentMopNumber,'Value',1,'fontsize',18,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'CurrentMopNumber=str2num(Vnav(2).String);PlotMobileLidarVolumes');

Vnav(3)=uicontrol(VolFig,'style','text','position',[200 690 200 25],...
    'string','VolumeElevCutoff','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',18);

Vnav(4)=uicontrol(VolFig,'style','edit','position',[410 690 100 25],...
    'string',num2str(VolumeElevCutoff),'Value',1,'fontsize',18,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'VolumeElevCutoff=str2num(Vnav(4).String);PlotMobileLidarVolumes');

Vnav(5)=uicontrol(VolFig,'style','text','position',[520 690 300 25],...
    'string','MinFractionSurveysWithData','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',18);

Vnav(6)=uicontrol(VolFig,'style','edit','position',[830 690 100 25],...
    'string',num2str(MinFractionSurveysWithData),'Value',1,'fontsize',18,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'MinFractionSurveysWithData=str2num(Vnav(6).String);PlotMobileLidarVolumes');

Vnav(7)=uicontrol(VolFig,'style','text','position',[940 690 350 25],...
    'string','MinFractionalCoverageAboveElevCutoff','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',18);

Vnav(8)=uicontrol(VolFig,'style','edit','position',[1300 690 100 25],...
    'string',num2str(MinFractionalCoverageAboveElevCutoff),'Value',1,'fontsize',18,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'MinFractionalCoverageAboveElevCutoff=str2num(Vnav(8).String);PlotMobileLidarVolumes');

PlotMobileLidarVolumes

