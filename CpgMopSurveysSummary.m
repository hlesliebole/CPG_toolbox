% CpgMopSurveysSummary.m
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
% "USACE" = USACE Shoals airborne topobathy Lidar (2009,2014)
% "CCC" = CA Coastal Conservancy funded 2009-2011 airborne Lidar
% "USGS" = USGS funded Lidar (2016)
% "KMair" = SIO Melville Group pre/post-El Nino Lidar (2015,2016)
% "iG8wheel" = SIO RTK GPS wheelie surveys
% "RTKdrone" = SIO RTK drone "structure from motion" surveys
%

fprintf('Loading struct array Survey from SurveyMasterListWithMops.mat...\n')
load('SurveyMasterListWithMops.mat','Survey');

fprintf('\nSurvey structural array:\n')
Survey

fprintf('\nFirst survey, Survey(1), meta data:\n')
Survey(1)

% number of surveys in the CPG MOP data set
NumSurveys=size(Survey,2);

% number of Mops with at least 1 survey 
NumMops=numel(unique([Survey.NearestMops]));

% number of unique survey dates 
NumDates=numel(unique([Survey.Datenum]));

fprintf('\nCPG MOPs includes %i surveys and %i unique Mop areas on %i unique dates.\n',...
    NumSurveys,NumMops,NumDates)

%% list the different kinds or sources of survey data
fprintf('\nSurvey Source    Number of Surveys   Number of Mops with data\n')

% unique survey sources
SurveySources=unique({Survey.Source},'stable'); 

for s=SurveySources
    % Survey indices matching source source
    idx=find(strcmp({Survey.Source},s)); 
    % number of mops with data from this survey source
    NumMops=numel(unique([Survey(idx).NearestMops])); 
    fprintf('%10s %15i %18i\n',string(s),numel(idx),NumMops)
end

%% example to list historical surveys associated with specific Mop number

% find all the Survey indices that have at least one data point that is 
%  nearest to the SIO Pier (Mop 513). 
ndx=find(cellfun(@(v)any(v(:) == 513),{Survey.NearestMops}));

% list all of them
fprintf('\nExample: SIO Pier Mop 513 Surveys:\n')
horzcat(string({Survey(ndx).Source}'),datestr([Survey(ndx).Datenum]'))