
%---------------------------------------------------
% updating daily iG8wheel subaerial transect surveys
%---------------------------------------------------

UpdateiG8wheel.m  % main program

% 1. edit the data file name and mop meta info in 
UpdateiG8wheelSAmatfiles.m
% and then run to update the SA files

% 2. update the SG files for the same mop range 
%  by editing the mop range in, and then run
AddLatestSGmatfiles.m

% 3. update the SM and GM files for the same mop range 
%  by editing the mop range in, and running
BuildSMmatfiles.m
%  this takes some time to run because it completely rebuilds
%  the GM files as well as the SM files.

% edit to select Mop number and number o recent profiles to plot.
PlotLatestSurveyProfiles

%---------------------------------------------------

% update latest iG8wheel survey
MakeiG8wheelSurveyList
% make struct array Survey

% check for struct Survey.File names in 
%  Survey.MopStart to Survey.MopEnd Mop SA matfiles 

% if not included yet, update SA, SG , and SM files
%  for GM array, update the month mean element and then 
%   use that to calculate the global mean from all the
%   month means.
%---------------------------------------------------


MakeSurveyMasterList.m

makes the struct array Survey 

Survey = 

  1×NsurveyFiles struct array with fields:

    DataSetNum
    File
    Datenum
    Source
    Type
    Year
    Month
    Day
    MopStart
    MopEnd


that contains info for all known processed survey files on reefbreak1.

% example code to find survey file names that are not included in the
%  present 1-m averaged survey SA struct array for a cpg mop.
idx=find(~ismember({Survey.File},{SA.File}) &...
    [Survey.MopStart] <= SA(1).Mopnum &...
    [Survey.MopEnd] >= SA(1).Mopnum);

%---------------------------------------

BuildCpgSurveyList.m
CompareMasterLists.m


BuildSAmatfilesMopRange.m

