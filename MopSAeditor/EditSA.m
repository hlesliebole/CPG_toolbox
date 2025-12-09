%
%  Main MATLAB script for editing SA files
%

% put subdirectory with the other mfiles in path

% clearvars
% close all


% initialize stuff
ns=0;
nb=0;
nn=0;

dview='3';

MainFig=figure('position',[250  103  1183 694]);

% default start is 654 fletcher cove if CurrentMopNumber does nor exist
if ~exist('CurrentMopNumber','var')
CurrentMopNumber=654;
end

if ~exist('CurrentSurveyNumber','var')
CurrentSurveyNumber=1;
end

% get the SA and QC struct arrays for the starting Mop
GetMopMobileLidarSA

