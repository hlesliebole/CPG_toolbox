%% CpgMopSurveysSummary.m
%
% The mat file SurveyMasterListWithMops.mat contains a struct array
% called "Survey" with the following meta data fields:
%
%     Survey.File     - full path to the processed survey file name on reefbreak
%     Survey.Bytes    - size of processed survey file
%     Survey.Datenum  - date of the actual survey
%     Survey.Source   - source or kind of survey (see list below)
%     Survey.Type     - a partial path survey file name
%     Survey.Year     - year of survey (for easier searching by year) 
%     Survey.Month    - month of survey (for easier searching by month)
%     Survey.Day      - day of survey  (for easier searching by day)
%     Survey.MopStart - starting Mop in file name's mop range
%     Survey.MopEnd   - ending Mop in file name's mop range
%     Survey.NearestMops  - list of actual Mop numbers that have data in the file
%     Survey.FileDatenum  - the creation date of the survey file 
%
% ----- Survey sources ------------
%
% "Gps" =  SIO in-situ Surveys (atv,dolly,jetski)
% "UTAir" = SIO UTexas Lidar (2000-2010) [including 1997 & 1998 NOAA/USGS ATM]
% "Trk" = SIO Truck Reigl Lidar (2018-)
% "TrkMR" = SIO Truck MiniRanger Lidar
% "AtvMR" = SIO ATV MiniRanger Lidar
% "DrnMR" = SIO Drone MiniRanger Lidar
% "UCACE" = UCACE Shoals airborne topobathy Lidar (2009,2014)
% "CCC" = CA Coastal Conservancy funded 2009-2011 airborne Lidar
% "USGS" = USGS funded Lidar (2016)
% "KMair" = SIO Melville Group pre/post-El Nino Lidar (2015,2016)
% "iG8wheel" = SIO RTK GPS wheelie surveys
% "RTKdrone" = SIO RTK drone "structure from motion" surveys
%

clearvars

CpgDefineMopPath

fprintf('Loading struct array Survey from SurveyMasterListWithMops.mat...\n')
load('SurveyMasterListWithMops.mat','Survey');

% find all the Mop numbers with survey data
MopNums=unique([Survey.NearestMops]);
MopNums=sort(MopNums);

% the CA files have to be updated first because the SG gridding uses
%  CA data from the surrounding Mops.

for m=654%MopNums
  fprintf(' ------ Mop %i -----\n',m); 

% first just read in the smaller CAfiles inventory struct array for this mop 
%  to check for changes
 
% read in the CA file if it exists
matfile=[ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ];
if exist(matfile,'file')
 fprintf('Loading CAfiles from %s\n',matfile);
 load([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CAfiles');
 NumCAsurveys=size(CAfiles,2);
 if size(CAfiles,2) == 0 % check if CAfiles is missing
     load([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA');
     CAfiles=struct('File',{CA.File},'FileDatenum',...
            num2cell([CA.FileDatenum]));
     save([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA','CAfiles');
 end
else
    fprintf('%s does not exist, making empty CAfiles struct array.\n',matfile);
    CAfiles(1).File=' ';
    CAfiles(1).FileDatenum=0;
end

% % find all the Survey indices that have at least one data point that is 
% %  nearest to the current Mop. 

ndx=find(cellfun(@(v)any(v(:) == m),{Survey.NearestMops}));
fprintf('The current Survey struct array has %i surveys for this Mop.\n',...
    numel(ndx))

% combine file datenums and file names into single strings for comparison
% of both the file names and creation dates at the CAme time to see if
% anything has changed

   Adf=strcat(string([CAfiles.FileDatenum]),{CAfiles.File}); % SG string
   Sdf=strcat(string([Survey(ndx).FileDatenum]),{Survey(ndx).File}); % Survey string
  
   ndel=find(~ismember(Adf,Sdf)); % CA file index that is not a member of Survey
   nadd=find(~ismember(Sdf,Adf)); % Survey file index that is not a member of CA
   
% ----------
% If any changes detected, load the full CA struct array and make changes
% ----------

if numel(ndel) > 0 || numel(nadd) > 0

 if exist(matfile,'file')
  load([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA');
  NumCAsurveys=size(CA,2);
 else
   CA=[];
   NumCAsurveys=0;
 end
    
% delete CA survey file entries that no longer exist or have changed on
% reefbreak
 if numel(ndel) > 0 && size(CA,2) > 0
  fprintf('Deleting %i CA surveys that have changed.\n',numel(ndel));
  CA(ndel)=[];
  % save the revised Mops CA file with old entries deleted
  CAfiles=struct('File',{CA.File},'FileDatenum',num2cell([CA.FileDatenum]));
  save([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA','CAfiles');
 end


pause

% add new or changed surveys to Mop CA files (this Mop and any other Mops
%  the survey covers.
 if numel(nadd) > 0
  fprintf('Adding %i CA surveys that are new or have changed.\n',numel(nadd));
  % loop through new survey files indices in Survey
  for n=nadd
    [filepath,name,ext] = fileparts(Survey(ndx(n)).File);
    fprintf('Adding survey %s : %s%s\n',Survey(ndx(n)).Source,name,ext)
    % update all the Mop CA files associated with this survey file
    % check to make sure group and drone folders on reefbreak are mounted
    %if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
    if exist('/volumes/group','dir') 
     % survey file update returns the struct arrCAy Survey that may have
     %  had its Survey.NearestMops fields modified depending on whether
     %  a Mop contained data after being reduced to 1m spatial averages.
     Survey=CpgUpdateCAmatfilesSingleSurvey(Survey,ndx(n));
    else
     %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
     fprintf('***** ERROR ******\n/volumes/group is not mounted.\nStopping.\n')
     break
    end
  end
  
 end

else
    fprintf('No changes detected.\n');
end 

end

% save Survey master list in case some Survey.NearestMops fields have been
% changed.
save('SurveyMasterListWithMops.mat','Survey');

% idx=[]; 
% for s=NumCAsurveys:-1:1
%     % find CAme file name and creation date in Survey if it exists
%     idx=find(strcmp({Survey(ndx).File},CAfiles(s).File) &...
%         [Survey(ndx).FileDatenum] == CAfiles(s).FileDatenum);
%     if numel(idx) == 0 
%         break
%     end
% end
% 
% if numel(idx) == 0 
%     
%  fprintf('Changes found. Loading CA struct array from %s\n',matfile);
%  load([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA');
%  NumCAsurveys=size(CA,2); 
%     
% % counter of suevey entries in CA that need to be deleted
% ndel=0; 
% 
% % logical array of new Survey indices for this mop
% newndx=ones(size(ndx)); % initialize as all 1 = new surveys
% 
% % Loop through survey file entries in CA.
% 
% %CAfiles=struct('File',{CA.File},'FileDatenum',num2cell([CA.FileDatenum]));
% 
% % Delete any CA entries with filenames and/or creation dates that 
% % don't match any entries in the master Survey struct array.
% 
% % Loop from last to first so entries can be delated without
% % changing indices of CA entries not checked yet
% for s=NumCAsurveys:-1:1
%     % find CAme file name and creation date in Survey if it exists
%     idx=find(strcmp({Survey(ndx).File},CA(s).File) &...
%         [Survey(ndx).FileDatenum] == CA(s).FileDatenum);
%     if numel(idx) == 0 
%         ndel=ndel+1;  % survey in CA no longer exists in Survey.
%         [filepath,name,ext] = fileparts(Survey(s).File);
%         fprintf('Deleting survey %s : %s%s\n',Survey(s).Source,name,ext)
%         CA(s)=[]; % delete the CA entry
%     else
%         % survey already exists in CA
%         newndx(idx)=0; % flag as 0 = not new  
%     end
% end
% 
% fprintf('Deleted %i of the %i CA surveys\n',ndel,NumCAsurveys)
% fprintf('Survey has %i new surveys for CA\n',numel(find(newndx == 1)))
% 
% % save the current Mops CA file with old entries deleted
% CAfiles=struct('File',{CA.File},'FileDatenum',num2cell([CA.FileDatenum]));
% save([ mpath 'M'  num2str( m , '%5.5i' )  'CA.mat' ],'CA','CAfiles');
% 
% % loop through the new Survey entries and add the survey data to both 
% %   this CA file and the CA files of any other mops that the survey 
% %   includes  (ie. only want to read and process the survey file once).
% 
% ndx=ndx(newndx == 1);
% 
% % loop through new survey files indices in Survey
% for n=ndx
%     [filepath,name,ext] = fileparts(Survey(n).File);
%     fprintf('Adding survey %s : %s%s\n',Survey(n).Source,name,ext)
%     % update all the Mop CA files associated with this survey file
%     % check to make sure group and drone folders on reefbreak are mounted
%     %if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
%     if exist('/volumes/group','dir') 
%      UpdateCAmatfilesSingleSurvey(Survey,n)
%     else
%      %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
%      fprintf('***** ERROR ******\n/volumes/group is not mounted.\nStopping.\n')
%      break
%     end
% end
% 
% else
%     fprintf('No changes found. \n');
% end
% 
%     %if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
%     if exist('/volumes/group','dir') 
%     
%     else
%      %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
%      fprintf('***** ERROR ******\n /volumes/group  is not mounted.\nStopping.\n')
%      break
%     end
% 
% end
% 
% % need to constantly check if connected to reefbreak and stop the
% %  update if not. The code is designed to be restarted/rerun without
% %  problems if a disconnect occurs.
% % 
% % 
% % fprintf('\nSurvey structural array:\n')
% % Survey
% % 
% % fprintf('\nFirst survey, Survey(1), meta data:\n')
% % Survey(1)
% % 
% % % number of surveys in the CPG MOP data set
% % NumSurveys=size(Survey,2);
% % 
% % % number of Mops with at least 1 survey 
% % NumMops=numel(unique([Survey.NearestMops]));
% % 
% % % number of unique survey dates 
% % NumDates=numel(unique([Survey.Datenum]));
% % 
% % fprintf('\nCPG MOPs includes %i surveys and %i unique Mop areas on %i unique dates.\n',...
% %     NumSurveys,NumMops,NumDates)
% % 
% % %% list the different kinds or sources of survey data
% % fprintf('\nSurvey Source    Number of Surveys   Number of Mops with data\n')
% % 
% % % unique survey sources
% % SurveySources=unique({Survey.Source},'stable'); 
% % 
% % for s=SurveySources
% %     % Survey indices matching source source
% %     idx=find(strcmp({Survey.Source},s)); 
% %     % number of mops with data from this survey source
% %     NumMops=numel(unique([Survey(idx).NearestMops])); 
% %     fprintf('%10s %15i %18i\n',string(s),numel(idx),NumMops)
% % end
% % 
% % %% example to list historical surveys associated with specific Mop number
% % 
% % % find all the Survey indices that have at least one data point that is 
% % %  nearest to the SIO Pier (Mop 513). 
% % ndx=find(cellfun(@(v)any(v(:) == 513),{Survey.NearestMops}));
% % 
% % % list all of them
% % fprintf('\nExample: SIO Pier Mop 513 Surveys:\n')
% % horzcat(string({Survey(ndx).Source}'),datestr([Survey(ndx).Datenum]'))