% GetAllMopsContourHistoryV2.m

% 1. Uses nearest profile point method to construct profiles from 1m avg
%    survey data. Settings Xtol and Ytol control data distance/gap
%    tolerances when making profiles.
%
% 2. Defines a truck lidar back beach point on each Mop transect
%
% 3. Removes points seaward of the profile minima (assumed to be swash 
%     or other errors)
%
% 4. Finds profile point/location (in Truck back beach coordimates) that 
%      is closest in elevation to a desired contour.  The point is considered
%      valid if its elev is within dztol of the contour elev.
%   
% Saves survey date, mop number, msl and mhw lat lon intercepts with the
%  mop transect, and msl and mhw distances in meters from the mop back
%  beach points, in the table array MopShoreline.

%  The MopShoreline table array is saved in MopShoreline.mat

%mpath='/volumes/group/MOPS/'; % file path to Mop files
mpath='/Users/William/Desktop/MOPS/'; % file path to Mop files
Ytol=50; % 50m alongcoast nearest point profile tolerance
Xtol=10; % 5m  nearest profile cross shore gap interpolation tolerance
dztol=0.1; % nearest profile point elev diff from contour elev tolerance

Mop1=3; % MX border 
Mop2=11594; % OR border

% path to GetNearestPointsProfiles.m function
eval(['addpath ' mpath 'toolbox'])

%%--------------------------------
%% Make Shoreline Data Table Array
%%--------------------------------

% create table with 6 columns of Data Set Number, (YYYYMMDD) SurveyDate,
% MopNum, MhwEast, MhwNorth, MhwDist
% where MhwDist is distance (m) from the mop back 
% beach point to mean high water (1.344m navd88).
% DataSetNum=[];
DateNum=[];
Year=[];
Month=[];
Day=[];
MopNum=[];
TrkX0=[];
% MhwEast=[];
% MhwNorth=[];
%MhhwDist=[];
MhwDist=[];
% MslEast=[];
% MslNorth=[];
MslDist=[];

% MopShoreline=table(MopNum,DataSetNum,DateNum,Year,Month,Day,...
%     MhwEast,MhwNorth,MhwDist,MslEast,MslNorth,MslDist);
MopShoreline=table(MopNum,TrkX0,DateNum,Year,Month,Day,MhwDist,MslDist);

% define tide elevations navd88
mllw=-0.058;mlw=0.218;msl=0.774;mhw=1.344;mhhw=1.566;hat=2.119;

%% loop through all possible mop numbers

 for MopNumber=Mop1:Mop2
  SA=[];
    matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SA.mat' ];
    if exist(matfile,'file')
     fprintf('Getting SA structat array for Mop %i of %i\n',MopNumber,Mop2)
     load(matfile,'SA');
    end
 

if size(SA,2) > 0 && ~isempty(SA(1).Mopnum) % check for empty struct array
fprintf('%i surveys.\n',size(SA,2))

%% (1) Turn 1m avg survey data into nearest point mop transect profiles
% X1D is a universal vector of 1m spatial res xshore profile locations 
% in the Mop back beach point coordinate frame.  
% Z1D(n,m) matrix of profile elevations with n=number of surveys in SA, 
% and m=size(X1D) is xshore profile coords.
[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);
% find the truck lidar surveys 
trk=find(strcmp({SA.Source},'Trk'));

if ~isempty(trk)
    xback=[];
% get the profile back beach location (in Mop transect x coords) for each 
% truck survey
 for n=1:numel(trk)
     xb=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
     if isempty(xb);xb=NaN;end
     xback(n)=xb;
 end

%% (2) Define the truck back beach (in Mop x coords) as the median truck 
%  back beach point on the mop transect
 X0BeachOnly=nanmedian(xback);
 fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)
 if isnan(X0BeachOnly)
     X0BeachOnly=0; % didn't find a truck back beach point, stay with Mop point
 else
     % truck back beach point defined:
     %  shift and reduce X1D locations and Z1D range 
     %  based on truck beach-only back beach location
     idx=find(X1D >= X0BeachOnly);    
     X1D=0:numel(idx)-1; % new xshore distances relative to X0BeachOnly
     Z1D=Z1D(:,idx); % keep elev data seaward of truck back beach
 end

else
    fprintf('No Truck Data, using Mop back beach point.\n',X0BeachOnly)
    X0BeachOnly=0; % stay with mop point
end
    
%% (3) Remove any elevations seaward of the min profile elevation
%  (lidar swash filter)
[zmin,imin]=min(Z1D');
for n=1:size(Z1D,1)
    if imin(n) > 0
        Z1D(n,imin(n)+1:end)=NaN; % set everything seaward of the minima to NaN
    end
end

%% (4) Find nearest point to MSL and MHW contours

[zmin,iminMSL]=min(abs(Z1D-msl)');
iminMSL(isnan(zmin))=NaN;
iminMSL(zmin > dztol)=NaN;
fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(iminMSL))),msl)

[zmin,iminMHW]=min(abs(Z1D-mhw)');
iminMHW(isnan(zmin))=NaN;
iminMHW(zmin > dztol)=NaN;
fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(iminMHW))),mhw)

%if numel(find(~isnan(imin))) > 0

% add row to shoreline table based on mops with Mhw contour location
   for i=1:size(SA,2)
     Shoreline.MopNum=MopNumber;
     Shoreline.TrkX0=X0BeachOnly;
     %Shoreline.DataSetNum=Survey(n).DataSetNum;
     Shoreline.DateNum=SA(i).Datenum;
     Shoreline.Year=year(datetime(SA(i).Datenum,'convertfrom','datenum'));
     Shoreline.Month=month(datetime(SA(i).Datenum,'convertfrom','datenum'));
     Shoreline.Day=day(datetime(SA(i).Datenum,'convertfrom','datenum')); 
%      Shoreline.MhwEast=cx(i);
%      Shoreline.MhwNorth=cy(i);
     if isnan(iminMHW(i))
         Shoreline.MhwDist=NaN;
     else
         Shoreline.MhwDist=X1D(iminMHW(i)); % distance relative to truck back beach
     end
     
     if isnan(iminMSL(i))
         Shoreline.MslDist=NaN;
     else
         Shoreline.MslDist=X1D(iminMSL(i)); % distance relative to truck back beach
     end
     
     MopShoreline = [MopShoreline;struct2table(Shoreline)];
   end    

end

 end
 
save MopShorelineV2.mat MopShoreline
% 
% 
% 
% 
% 
% %    *Uses m-script mopcontour.m
% % 
% 
% % reads in all avialble historical beach survey data and derives estimates 
% % of the mhw mop contour locations at each mop transect, and its 
% % distance from each mop back beach point.
% 
% % requires mounting the reefbreak1 /project/group volume first, using
% % Finder -> Go -> Connect to Server 
% % The reefbreak1 group data will then appear under /Volumes/group/ 
% % on your mac.
% 
% % reads both atv-jumbo llzts.navd88 ascii files and truck lidar tif
% % 1-m resolution survey data files
% 
% % Saves survey date, mop number, msl and mhw lat lon intercepts with the
% %  mop transect, and msl and mhw distances in meters from the mop back
% %  beach points, in the table array MopShoreline.
% 
% %  The MopShoreline table array is saved in MopShoreline.mat
% 
% 
% % get definition of survey data sets to include
% Step1defineDataSets
% 
% % get master list of survey files
% %Step2getSurveyFileList
% addpath /Users/William/Desktop/MOPS/toolbox
% load SurveyMasterListWithMops.mat  % quick option to load prexisting survey struct array
% 
% %-------------------------------
% % Load Mop information
% %-------------------------------
% % path to MOP transect definition text file 
% mopurl='http://cdip.ucsd.edu/MOP_v1.1/CA_v1.1_transect_definitions.txt';
% 
% fprintf('Loading Mop Transect Location Info from:\n %s \n',mopurl)
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
% 
% %--------------------------------
% % Make Shoreline Data Table Array
% %--------------------------------
% 
% % create table with 6 columns of Data Set Number, (YYYYMMDD) SurveyDate,
% % MopNum, MhwEast, MhwNorth, MhwDist
% % where MhwDist is distance (m) from the mop back 
% % beach point to mean high water (1.344m navd88).
% DataSetNum=[];
% DateNum=[];
% Year=[];
% Month=[];
% Day=[];
% MopNum=[];
% MhwEast=[];
% MhwNorth=[];
% MhwDist=[];
% MslEast=[];
% MslNorth=[];
% MslDist=[];
% 
% MopShoreline=table(MopNum,DataSetNum,DateNum,Year,Month,Day,...
%     MhwEast,MhwNorth,MhwDist,MslEast,MslNorth,MslDist);
% 
% % define tide elevations navd88
% mllw=-0.058;mlw=0.218;msl=0.774;mhw=1.344;mhhw=1.566;hat=2.119;
% 
% %------------------------------------------
% % Loop thru survey master list struct array
% %   and add each mop mean mhw contour distance to table
% %------------------------------------------
% 
% j=0; % initializeoutput table row counter
% 
% for n=1:size(Survey,2)
% %for n=1:2
% fprintf('Survey number: %g \n',n)
%  
%  % read in x,y,z data in UTM coords
%  [xutm,yutm,z]=readSurveyFileUTM(Survey(n).source,Survey(n).file); 
%  
%  % skip if all the survey day is above or all below mhw
%  if( mhw > min(z) && mhw < max(z) )
% 
%  % get mop mhw contour info
%    [cmop,cx,cy,cdist]=mopcontourUTM(Mop,xutm,yutm,z,mhw);
%    
%    
%    % get mop msl contour info
%    cmop2=[];
%     % skip if all the survey day is above or all below msl
%     if( msl > min(z) && msl < max(z) )
%      [cmop2,cx2,cy2,cdist2]=mopcontourUTM(Mop,xutm,yutm,z,msl);
%     end
%  
%  % add row to shoreline table based on mops with Mhw contour location
%    for i=1:length(cmop)
%      Shoreline.MopNum=cmop(i);
%      Shoreline.DataSetNum=Survey(n).DataSetNum;
%      Shoreline.DateNum=Survey(n).datenum;
%      Shoreline.Year=Survey(n).year;
%      Shoreline.Month=Survey(n).month;
%      Shoreline.Day=Survey(n).day; 
%      Shoreline.MhwEast=cx(i);
%      Shoreline.MhwNorth=cy(i);
%      Shoreline.MhwDist=cdist(i);
%      % add msl location if one was found
%      imsl=find(cmop2 == cmop(i));
%      if ~isempty(imsl)
%          Shoreline.MslEast=cx2(imsl);
%          Shoreline.MslNorth=cy2(imsl);
%          Shoreline.MslDist=cdist2(imsl);
%      else % otherwise fill msl info with NaNs
%          Shoreline.MslEast=NaN;
%          Shoreline.MslNorth=NaN;
%          Shoreline.MslDist=NaN;
%      end
%          
%      MopShoreline = [MopShoreline;struct2table(Shoreline)];
%    end
%  
%  end
%  
% end
% 
% % sort table by datenum then by mop number
% 
%  MopShoreline = sortrows(MopShoreline,3);
%  MopShoreline = sortrows(MopShoreline,1);
% 
% fprintf(1,' Saving MopShoreline table array to MopShoreline.mat \n')
% 
% save MopShoreline.mat MopShoreline
% 
