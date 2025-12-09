% Makes "DataSets" struct array with key data set information

% After you run this script, the wildcard file search path for each
% source of data will be in the struct variable element
%
%    DataSets(n).filepath
%
% where n is one of the 16 sources listed below.



% A survey type can have more than one DataSet definition depending on file
%  naming conventions.

%    DataSets(n).name : Data set name (Gps,UTAir,Truck)
%
%    DataSets(n).title : Character string description of data set
%
%    DataSets(n).filehead : Defines the portion of the filepath, usually up to 
%    encountering the YYYYMMDD_MOPSTART_MOPEND portion of the file folder 
%    or filename string. This could be either a folder name or a folder name 
%    plus a portion of the file name depending on the data source folder and
%    file name conventions.  This is handy for parsing out the date and mop
%    range associated with a survey data file.
%
%    DataSets(n).filepath :  the actual wildcard search string to use
%    when looking for files with the matlab dir command
%

% Settings for mounted reefbreak1 data directories
% Check if this is a pc to set the reefbreak1 group directory 

if ispc
    group='Y:\';  % Adam's PC mount point
else
    group='/Volumes/group/'; % Mac mount point
    drone='/Volumes/drone/'; % Mac path for
end

%-------------------------------------------------------------------------
% 1. SIO Truck LiDAR (VMZ2000 system)
%
% group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkFiles='/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(1).name='Trk';
DataSets(1).title='SIO Truck Lidar Surveys';
DataSets(1).filehead=...
fullfile(group,'LiDAR','VMZ2000_Truck','LiDAR_Processed_Level2',filesep);
DataSets(1).filepath=...
fullfile(group,'LiDAR','VMZ2000_Truck','LiDAR_Processed_Level2',...
'*','Beach_And_Backshore','2*ground.tif');

%-------------------------------------------------------------------------
% 2. SIO Truck MiniRanger 
%
% group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkMRFiles='/Volumes/group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(2).name='TrkMR';
DataSets(2).title='SIO Truck MiniRanger Surveys';
DataSets(2).filehead=...
fullfile(group,'LiDAR','MiniRanger_Truck','LiDAR_Processed_Level2',filesep);
DataSets(2).filepath=...
fullfile(group,'LiDAR','MiniRanger_Truck','LiDAR_Processed_Level2',...
'*','Beach_And_Backshore','2*ground.tif');


%-------------------------------------------------------------------------
% 3. SIO ATV MiniRanger 
%
% group/LiDAR/MiniRanger_ATV/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% AtvMRFiles='/Volumes/group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(3).name='AtvMR';
DataSets(3).title='SIO ATV MiniRanger Surveys';
DataSets(3).filehead=...
fullfile(group,'LiDAR','MiniRanger_ATV','LiDAR_Processed_Level2',filesep);
DataSets(3).filepath=...
fullfile(group,'LiDAR','MiniRanger_ATV','LiDAR_Processed_Level2',...
'*','Beach_And_Backshore','2*ground.tif');


%-------------------------------------------------------------------------
% 4. SIO Truck LiDAR (VMQLZ system)
%
% group/LiDAR/VMQLZ_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkFiles='/Volumes/group/LiDAR/VMQLZ_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(4).name='Trk';
DataSets(4).title='SIO Truck Lidar Surveys';
DataSets(4).filehead=...
fullfile(group,'LiDAR','VMQLZ_Truck','LiDAR_Processed_Level2',filesep);
DataSets(4).filepath=...
fullfile(group,'LiDAR','VMQLZ_Truck','LiDAR_Processed_Level2',...
'*','Beach_And_Backshore','2*ground.tif');

