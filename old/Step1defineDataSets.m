% Step1defineDataSets.m

% Makes "DataSets" struct array with key data set information needed
%  to fine the data files.

%
% This version defines 3 San Diego County topobathy survey types
%  that are currently stored in reefbreak1/group
%
% - SIO GPS Surveys (atv,dolly,jetski)
% - SIO UTexas Lidar (2000-2010) [including 1997 & 1998 NOAA/USGS ATM]
% - SIO Truck Lidar (2018-)

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

% Seetings for mounted reefbreak1 data directories
% Check if this is a pc to set the reefbreak1 group directory 

if ispc
    group='Y:\';  % Adam's PC mount point
else
    group='/Volumes/group/'; % Mac mount point
end

%-------------------------------------------------------------------------
% 1. atv-jumbo (Gps) survey data 
%

% Note 1: different file naming conventions requires 2 data set definitions/searches 
% to find all the Gps survey files
%
% Prior to 20040206, there is just a *.adj2011 file.
%
% Starting with  topobathy/20040206_00520_00597_posncex_atv_alongshore
% The folder contain both a  *.adj2011 file and *.navd88 file. file
% processing logic is needed to choose the navd88 file when both
% have been found with the separate searches.
% 
% Starting with topobathy/20111010_00520_00598_posncex_jumbo
% The folders only have a *.navd88 file
%
% Note 2: The .llz.navd88 files can have 3 (lat,lon,z), 4 (lat,lon,z,time), 
% 5 (lat,lon,z,time,substrate_code), or 
% 7 (lat,lon,Nutm,Eutm,z,time,substrate_code) columns. Not an issue when
% defining the data set paths here, but requires flexible code when reading
% the data files themselves as there is no way to tell someof these apart
% from the file name itself.

DataSets(1).name='Gps';
DataSets(1).title='SIO In Situ GPS Surveys';
%DataSets(1).filehead='/Volumes/group/topobathy/';
DataSets(1).filehead=fullfile(group,'topobathy',filesep);
%DataSets(1).filepath='/Volumes/group/topobathy/*/filtered*.ll*.adj2011';
DataSets(1).filepath=fullfile(group,'topobathy','*','filtered*.ll*.adj2011');

DataSets(2).name='Gps';
DataSets(2).title='SIO In Situ GPS Surveys';
%DataSets(2).filehead='/Volumes/group/topobathy/';
DataSets(2).filehead=fullfile(group,'topobathy',filesep);
%DataSets(2).filepath='/Volumes/group/topobathy/*/filtered*.ll*.navd88';
DataSets(2).filepath=fullfile(group,'topobathy','*','filtered*.ll*.navd88');

%-------------------------------------------------------------------------
% 2. SIO/UTexas Airborne LiDAR 
%
% group/LiDAR/LiDAR_Airborne/YYYYMMDD_MopStart_MopEnd_NoWaves_SDAirborne/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SDAirborne_beach_ground.txt
%

% AirFiles='/Volumes/group/LiDAR/LiDAR_Airborne/*/Beach_Only/*ground.txt';

% 3 columns  (Lon,Lat,z)

DataSets(3).name='UTAir';
DataSets(3).title='SIO-UT Airborne Lidar Surveys';
%DataSets(3).filehead='/Volumes/group/LiDAR/LiDAR_Airborne/';
DataSets(3).filehead=fullfile(group,'LiDAR','LiDAR_Airborne',filesep);
%DataSets(3).filepath='/Volumes/group/LiDAR/LiDAR_Airborne/*/Beach_Only/*ground.txt';
DataSets(3).filepath=...
    fullfile(group,'LiDAR','LiDAR_Airborne','*','Beach_Only','*ground.txt');


%-------------------------------------------------------------------------
% 3. SIO Truck LiDAR 
%
% group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkFiles='/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(4).name='Trk';
DataSets(4).title='SIO Truck Lidar Surveys';
%DataSets(4).filehead='/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/';
DataSets(4).filehead=...
fullfile(group,'LiDAR','VMZ2000_Truck','LiDAR_Processed_Level2',filesep);
%DataSets(4).filepath='/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/*/Beach_Only/2*ground.tif';
DataSets(4).filepath=...
fullfile(group,'LiDAR','VMZ2000_Truck','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');

%-------------------------------------------------------------------------
% 4. SIO Truck MiniRanger 
%
% group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkMRFiles='/Volumes/group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(5).name='TrkMR';
DataSets(5).title='SIO Truck MiniRanger Surveys';
DataSets(5).filehead=...
fullfile(group,'LiDAR','MiniRanger_Truck','LiDAR_Processed_Level2',filesep);
DataSets(5).filepath=...
fullfile(group,'LiDAR','MiniRanger_Truck','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');


%-------------------------------------------------------------------------
% 5. SIO ATV MiniRanger 
%
% group/LiDAR/MiniRanger_ATV/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% AtvMRFiles='/Volumes/group/LiDAR/MiniRanger_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(6).name='AtvMR';
DataSets(6).title='SIO ATV MiniRanger Surveys';
DataSets(6).filehead=...
fullfile(group,'LiDAR','MiniRanger_ATV','LiDAR_Processed_Level2',filesep);
DataSets(6).filepath=...
fullfile(group,'LiDAR','MiniRanger_ATV','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');


%-------------------------------------------------------------------------
% 6. SIO Drone MiniRanger 
%
% group/LiDAR/MiniRanger_Drone/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% DrnMRFiles='/Volumes/group/LiDAR/MiniRanger_Drone/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(7).name='DrnMR';
DataSets(7).title='SIO Drone MiniRanger Surveys';
DataSets(7).filehead=...
fullfile(group,'LiDAR','MiniRanger_Drone','LiDAR_Processed_Level2',filesep);
DataSets(7).filepath=...
fullfile(group,'LiDAR','MiniRanger_Drone','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');


DataSets