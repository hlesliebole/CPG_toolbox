% Makes "DataSets" struct array with key data set information

% After you run this script, the wildcard file search path for each
% source of data will be in the struct variable element
%
%    DataSets(n).filepath
%
% where n is one of the 16 sources listed below.

%
% San Diego County topobathy datasets
%
% 1.-2. SIO GPS Surveys (atv,dolly,jetski)
% 3. SIO UTexas Lidar (2000-2010) [including 1997 & 1998 NOAA/USGS ATM]
% 4. SIO Truck Lidar (VMZ2000 2017-7/2024)
% 5. SIO Truck MiniRanger 
% 6. SIO ATV MiniRanger 
% 7. SIO Drone MiniRanger 
% 8. USACE Shoals topobathy Lidar (2009,2014)
% 9. Coastal Conservancy 2009-2011 Lidar
% 10. USGS Lidar (2016)
% 11. SIO Melville EL Nino Lidar (2015,2016)
% 12. SIO RTK wheel Survey
% 13. SIO OCS wheel Survey (Oceanside Citizen Science)
% 14. SIO Drone Structure from Motion
% 15. USACE Shoals Topbathy Airborne LiDAR in CCD_Airborne Files
% 16. USGS Topo Airborne LiDAR in CCD_Airborne Files
% 17. Norbit Multibeam
% 18. SIO ScoutUltra Backpack
% 19. 1997-98 NASA (UTAir) Airborne LiDAR from Digitial Coast
% 20. SIO Truck LiDAR (VMQLZ 7/2024-)

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
    if exist('/project/group/','dir')
        group='/project/group/'; % iwa mount point
        drone='/project/drone/'; 
   
    else
        group='/Volumes/group/'; % Mac mount point
        drone='/Volumes/drone/'; % Mac path for
    end
end

%-------------------------------------------------------------------------
% 1.-2. atv-jumbo (Gps) survey data 
%

% Note 1: different file naming conventions requires 2 data set definitions/searches 
% to find all the Gps survey files

%
% Note 2: The .llz.navd88 files can have 3 (lat,lon,z), 4 (lat,lon,z,time), 
% 5 (lat,lon,z,time,substrate_code), or 
% 7 (lat,lon,Nutm,Eutm,z,time,substrate_code) columns. Not an issue when
% defining the data set paths here, but requires flexible code when reading
% the data files themselves as there is no way to tell from the file name
% itself.

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
% 3. SIO/UTexas Airborne LiDAR 
%
% group/LiDAR/LiDAR_Airborne/YYYYMMDD_MopStart_MopEnd_NoWaves_SDAirborne/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SDAirborne_beach_ground.txt
%

% AirFiles='/Volumes/group/LiDAR/LiDAR_Airborne/*/Beach_Only/*ground.txt';

% 3 columns  (Lon,Lat,z)

DataSets(3).name='UTAir';
DataSets(3).title='SIO-UT Airborne Lidar Surveys';
%DataSets(3).filehead='/Volumes/group/LiDAR/LiDAR_Airborne/';
DataSets(3).filehead=fullfile(group,'LiDAR','LiDAR_Airborne','Topo','Level_2',filesep);
%DataSets(3).filepath='/Volumes/group/LiDAR/LiDAR_Airborne/*/Beach_Only/*ground.txt';
% DataSets(3).filepath=...
%     fullfile(group,'LiDAR','LiDAR_Airborne','Level_2','*','Beach_Only','*ground.txt');
DataSets(3).filepath=...
    fullfile(group,'LiDAR','LiDAR_Airborne','Topo','Level_2','*','Beach_Only','*ground.txt');


%-------------------------------------------------------------------------
% 4. SIO Truck LiDAR (VMZ2000 system)
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
% 5. SIO Truck MiniRanger 
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
% 6. SIO ATV MiniRanger 
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
% 7. SIO Drone MiniRanger 
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

%-------------------------------------------------------------------------
% 8. USACE Shoals Topbathy Airborne LiDAR 
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_USACEShoals.llz
%
% 3 columns  (Lon,Lat,z)

DataSets(8).name='USACE';
DataSets(8).title='USACE Shoals Topobathy Lidar Surveys';
DataSets(8).filehead=...
fullfile(group,'BOR','DigitalCoastLidar',filesep);
DataSets(8).filepath=...
fullfile(group,'BOR','DigitalCoastLidar','2*USACEShoals.llz');

%-------------------------------------------------------------------------
% 9. CA Coastal Conservancy Airborne LiDAR (Fugro) 
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_CCClidar.llz
%
% 3 columns  (Lon,Lat,z)

DataSets(9).name='CCC';
DataSets(9).title='CA Coastal Conservancy Lidar Surveys';
DataSets(9).filehead=...
fullfile(group,'BOR','DigitalCoastLidar',filesep);
DataSets(9).filepath=...
fullfile(group,'BOR','DigitalCoastLidar','2*CCClidar.llz');

%-------------------------------------------------------------------------
% 10. USGS Airborne LiDAR 
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_USGSlidar.llz
%
% 3 columns  (Lon,Lat,z)

DataSets(10).name='USGS';
DataSets(10).title='USGS Lidar Surveys';
DataSets(10).filehead=...
fullfile(group,'BOR','DigitalCoastLidar',filesep);
DataSets(10).filepath=...
fullfile(group,'BOR','DigitalCoastLidar','2*USGSlidar.llz');

%-------------------------------------------------------------------------
% 11. Melville Airborne Lidar 
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_KMAir.llz
%
% 3 columns  (Lon,Lat,z)

DataSets(11).name='KMair';
DataSets(11).title='SIO Melville Group Lidar Surveys';
DataSets(11).filehead=...
fullfile(group,'BOR','DigitalCoastLidar',filesep);
DataSets(11).filepath=...
fullfile(group,'BOR','DigitalCoastLidar','2*KMair.llz');

%-------------------------------------------------------------------------
% 12. RTK iG8 wheel survey 
%
% /Volumes/group/ig8_wheel/YYYYMMDD_MopStart_MopEnd_*_wheel/filtered_cleanYYYYMMDD.llnezts.navd88
%
% 6 columns  (Lon,Lat,Eutm,Nutm,z,t)

DataSets(12).name='iG8wheel';
DataSets(12).title='RTK iG8 wheel Surveys';
DataSets(12).filehead=...
fullfile(group,'ig8_wheel',filesep);
DataSets(12).filepath=...
fullfile(group,'ig8_wheel','2*wheel','filtered*.navd88');

%-------------------------------------------------------------------------
% 13. Oceanside citizen wheel surveys 
%
% /Volumes/group/Oceanside/YYYYMMDD_MopStart_MopEnd_oceanside_wheel/filtered*.navd88
%
%  need to convert txt files to llz files in folders with
%
% cat X.txt | cut -f2-3,8 -d' ' > X.llz

DataSets(13).name='OCSwheel';
DataSets(13).title='Oceanside Citizen Surveys';
DataSets(13).filehead=...
fullfile(group,'Oceanside',filesep);
DataSets(13).filepath=...
fullfile(group,'Oceanside','2*oceanside_wheel*','filtered*.navd88');

%-------------------------------------------------------------------------
% 14. Drone Structure from Motion surveys 
%
% /Volumes/group/SfM/SfM_Processed_Level2/YYYYMMDD*/Beach_Only/YYYYMMDD_MopStart_MopEnd*.tif
%
% tif data file

DataSets(14).name='SfMdrone';
DataSets(14).title='SfM Drone Surveys';
DataSets(14).filehead=...
fullfile(drone,'SfM','torrey','2*',filesep);
DataSets(14).filehead=...
fullfile(group,'SfM','SfM_Processed_Level2',filesep);
DataSets(14).filepath=...
fullfile(group,'SfM','SfM_Processed_Level2',...
'*','Beach_Only','2*ground.tif');
%-------------------------------------------------------------------------
% 15. USACE Shoals Topbathy Airborne LiDAR in CCD_Airborne Files
%
% '/Volumes/group/LiDAR/CCD_Airborne/YYYYMMDD_MopStart_MopEnd_*USACE/Beach_Only/*.tif
%
DataSets(15).name='USACE';
DataSets(15).title='USACE Shoals Topobathy Lidar Surveys';
DataSets(15).filehead=...
fullfile(group,'LiDAR','CCD_Airborne',filesep);
DataSets(15).filepath=...
fullfile(group,'LiDAR','CCD_Airborne',...
'*USACE','Beach_Only','*.tif');
%-------------------------------------------------------------------------
% 16. USGS Topo Airborne LiDAR in CCD_Airborne Files
%
% '/Volumes/group/LiDAR/CCD_Airborne/YYYYMMDD_MopStart_MopEnd_*USGS/Beach_Only/*.tif
%
DataSets(16).name='USGS';
DataSets(16).title='USGS Topo Lidar Surveys';
DataSets(16).filehead=...
fullfile(group,'LiDAR','CCD_Airborne',filesep);
DataSets(16).filepath=...
fullfile(group,'LiDAR','CCD_Airborne',...
'*USGS','Beach_Only','*.tif');
%-------------------------------------------------------------------------
% 17. Norbit Multibeam
%
% '/Volumes/group/Multibeam/YYYYMMDD_MopStart_MopEnd_*/*multibeam.las
%
DataSets(17).name='Multibeam';
DataSets(17).title='Norbit Multibeam';
DataSets(17).filehead=...
fullfile(group,'Multibeam',filesep);
DataSets(17).filepath=...
fullfile(group,'Multibeam','2*multibeam','2*multibeam.las');

%-------------------------------------------------------------------------
% 18. SIO ScoutUltra Backpack
%
% group/LiDAR/ScoutUltra_Backpack/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% ScoutFiles='/Volumes/group/LiDAR/ScoutUltra_Backpack/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(18).name='Scout';
DataSets(18).title='SIO ScoutUltra Backpack Surveys';
DataSets(18).filehead=...
fullfile(group,'LiDAR','ScoutUltra_Backpack','LiDAR_Processed_Level2',filesep);
DataSets(18).filepath=...
fullfile(group,'LiDAR','ScoutUltra_Backpack','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');

% %-------------------------------------------------------------------------
% % 13. Torrey RTK Drone surveys 
% %
% % /Volumes/drone/data/torrey/YYYYMMDD/YYYYMMDD_MopStart_MopEnd*RTKdrone*.las
% %
% % las data file
% 
% DataSets(13).name='RTKdrone';
% DataSets(13).title='RTK Drone Surveys';
% DataSets(13).filehead=...
% fullfile(drone,'data','torrey','2*',filesep);
% %fullfile(drone,'data','torrey','20220216',filesep);
% DataSets(13).filepath=...
%     fullfile(drone,'data','torrey','2*','2*.las');

%-------------------------------------------------------------------------
% 19. USACE Shoals Topbathy Airborne LiDAR from Digitial Coast
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_USACEShoals.llz
%
% 3 columns  (Lon,Lat,z)
% 
% DataSets(19).name='USACE';
% DataSets(19).title='USACE Shoals Topobathy Lidar Surveys';
% DataSets(19).filehead=...
% fullfile(group,'BOR','DigitalCoastLidar',filesep);
% DataSets(19).filepath=...
% fullfile(group,'BOR','DigitalCoastLidar','2*USACEShoals.llz');

%-------------------------------------------------------------------------
% 19. 1997-98 NASA (UTAir) Airborne LiDAR from Digitial Coast
%
% /Volumes/group/BOR/DigitalCoastLidar/YYYYMMDD_MopStart_MopEnd_UTAir.llz
%
% 3 columns  (Lon,Lat,z)

DataSets(19).name='UTAir';
DataSets(19).title='SIO-UT Airborne Lidar Surveys';
DataSets(19).filehead=...
fullfile(group,'BOR','DigitalCoastLidar',filesep);
DataSets(19).filepath=...
fullfile(group,'BOR','DigitalCoastLidar','1*UTAir.llz');

%-------------------------------------------------------------------------
% 20. SIO Truck LiDAR (VMQLZ system)
%
% group/LiDAR/VMQLZ_Truck/LiDAR_Processed_Level2/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType/Beach_Only/YYYYMMDD_MopStart_MopEnd_NoWaves_SurveyType_ground.tif
%
%
% TrkFiles='/Volumes/group/LiDAR/VMQLZ_Truck/LiDAR_Processed_Level2/*/Beach_Only/*ground.tif';

DataSets(20).name='Trk';
DataSets(20).title='SIO Truck Lidar Surveys';
DataSets(20).filehead=...
fullfile(group,'LiDAR','VMQLZ_Truck','LiDAR_Processed_Level2',filesep);
DataSets(20).filepath=...
fullfile(group,'LiDAR','VMQLZ_Truck','LiDAR_Processed_Level2',...
'*','Beach_Only','2*ground.tif');

