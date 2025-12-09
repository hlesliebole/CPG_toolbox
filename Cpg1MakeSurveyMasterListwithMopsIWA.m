% addpath /volumes/group/mops/toolbox
% addpath /volumes/group/mops
clearvars
%CpgDefineMopPath

fprintf('\nAdding Nearest Mop Lists to Survey struct array.\n')

% check to make sure group and drone folders on reefbreak are mounted
%if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
if exist('/project/group','dir') 


load MopTableUTM.mat

% make standard list of current survey files on reefbreak1
MakeSurveyMasterList
% load that struct array as "CurrentSurvey" struct aray
load('SurveyMasterList.mat','Survey');
CurrentSurvey=Survey;
CurrentSurvey(1).NearestMops=[]; % add a nearest mops field to the struct array

% now load the last saved struct array of survey files, with their
%  associated .NearestMops with data lists, to struct array "Survey"
load('SurveyMasterListWithMops.mat','Survey');

% remove any past surveys that have no nearesr mop data
nsurv=size(Survey,2);
for n=nsurv:-1:1
    if numel(Survey(n).NearestMops) == 0
        fprintf('Removing %s\n from Survey for having no NearestMops values\n',Survey(n).File')
        Survey(n)=[];
    end
end

% Survey(1).FileDatenum=[];
% for n=1:size(Survey,2)
%     idx=find(strcmpi({CurrentSurvey.File},Survey(n).File));
%     if ~isempty(idx)
%       Survey(n).FileDatenum=CurrentSurvey(idx).FileDatenum;
%     end
% end

% find any survey file names in the previous Survey struct array
%  that are *not* in the CurrentSurvey struct array

idx1=find(~ismember(lower({Survey.File}),lower({CurrentSurvey.File})));

fprintf('-----------------------------------------------------\n')
fprintf('There are %i filenames in the previous Survey list that are not in CurrentSurvey.\n',numel(idx1))
fprintf('-----------------------------------------------------\n')

if numel(idx1) > 0
    for n=idx1
     fprintf('%s\n',Survey(n).File)
    end
end

% Now the opposite, find filenames in the CurrentSurvey struct array that
%  are not found in the Survey struct array 

idx2=find(~ismember(lower({CurrentSurvey.File}),lower({Survey.File})));

fprintf('-----------------------------------------------------\n')
fprintf('There are %i new filenames in CurrentSurvey that are not in Survey.\n',numel(idx2))
fprintf('-----------------------------------------------------\n')

if numel(idx2) > 0
    for n=idx2
     fprintf('%s\n',CurrentSurvey(n).File)
    end
end
fprintf('-----------------------------------------------------\n')

% see if the missing previous survey file are actually a renamed new file name
% for the same  date, source, and with the same number of bytes.  If so, change the 
%  Survey file name entry.

if numel(idx1) > 0
    for n=idx1
     fprintf(' Comparing: %s\n',Survey(n).File)
    if numel(idx2) > 0
    for n2=idx2
     fprintf('        To: %s\n',CurrentSurvey(n2).File)
%      fprintf('%i %i %i %i\n',...
%          abs(Survey(n).Datenum-CurrentSurvey(n2).Datenum),...
%          strcmp(Survey(n).Source,CurrentSurvey(n2).Source),...
%          Survey(n).Bytes,CurrentSurvey(n2).Bytes)
         
     if(abs(Survey(n).Datenum-CurrentSurvey(n2).Datenum) <= 1 ...
             && strcmp(Survey(n).Source,CurrentSurvey(n2).Source) == 1 ...
             && Survey(n).Bytes == CurrentSurvey(n2).Bytes)
         fprintf('          Data Match, Updating Filename in Survey %s\n',Survey(n).File)
            Survey(n).File=CurrentSurvey(n2).File;
     else
         %fprintf('          No Match.\n')
     end
    end
   end
    end
end

%-------------

% Now delete any remaining Survey entries that are not in CurrentSurvey
idx1=find(~ismember(lower({Survey.File}),lower({CurrentSurvey.File})));

if numel(idx1) > 0
    fprintf('Removing the following entries from Survey:\n')
    for n=fliplr(idx1)  % go from highest to lowest index
     fprintf('  %s\n',Survey(n).File)
     Survey(n)=[];
    end
end

% % Add new surveys and their nearest mops to Survey 

idx2=find(~ismember(lower({CurrentSurvey.File}),lower({Survey.File})));

if numel(idx2) > 0
    fprintf('Adding the following surveys to Survey:\n')
    nn=size(Survey,2);
    for n=idx2
      fprintf('%s\n',CurrentSurvey(n).File)
      nn=nn+1;
      Survey(nn)=CurrentSurvey(n);
     
      %[x,y,z,c,utmzone]=readSurveyFileUTM2(CurrentSurvey(n).Source,CurrentSurvey(n).File);

      [x,y,z,c]=readSurveyFileUTM(CurrentSurvey(n).Source,CurrentSurvey(n).File);
   
      if ~isempty(x) && numel(x) > 1
        [nmop]=XY2MopNumsV2(x,y,Mop);
          nmop=unique(sort(nmop));
      else
          nmop=[];
          fprintf('Warning: no NearestMops.\n')
      end
   
      Survey(nn).NearestMops=nmop;
      fprintf('Saving struct array Survey to SurveyMasterListWithMops.mat\n')
      save('SurveyMasterListWithMops.mat','Survey');

    end
end

% sort Survey by date before saving
T=struct2table(Survey);
sortedT = sortrows(T, 'Datenum');
Survey=table2struct(sortedT)';

fprintf('Saving struct array Survey to SurveyMasterListWithMops.mat\n')
save('SurveyMasterListWithMops.mat','Survey');

else
    %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
    fprintf('***** ERROR ******\nEither /volumes/group is not mounted.\nStopping.\n')
end

function  N=XY2MopNumsV2(xutm,yutm,Mop)

% calls DefineMopTransectUTM.m

% Returns a vector of the nearest mop transect numbers, N, based on the mop transect
% definition table, Mop, for the input x,y data point vectors.

%------------------------------------------------------------
%fprintf('Making approximate Mop match to start...\n')
% To narrow down the possible range of Mops in the survey area but avoid 
% excluding Mops at theedges of the survey region, set a distance tolerance, 
% tol, in meters around the min-max box of survey points.
tol=1000; % use 1000m to be safe
minx=min(xutm)-tol;maxx=max(xutm)+tol;
miny=min(yutm)-tol;maxy=max(yutm)+tol;
% find the Mops whose back beach or offshore point falls within
% the survey min-max + tolerance box.
mopidx=find( ((Mop.BackXutm >= minx & Mop.BackXutm <= maxx) | ...
        (Mop.OffXutm >= minx &  Mop.OffXutm <= maxx)) & ...
        ((Mop.BackYutm >= miny & Mop.BackYutm <= maxy) | ...
        (Mop.OffYutm >= miny &  Mop.OffYutm <= maxy)));
% reduce to possible start and end Mop numbers
MopStart=min(mopidx);MopEnd=max(mopidx);
% As first approximate matching of points to Mop transects, find
% closest of thse Mop transect back beach points (not the Mop transect line) 
% to each survey data point.
mtx=[];
mty=[];
mnum=[];
mtx=Mop.BackXutm(MopStart:MopEnd)'; 
mty=Mop.BackYutm(MopStart:MopEnd)'; 
mopnum=MopStart:MopEnd;
% match survey utm points to the nearest back beach points
[dp,NearIdx]=pdist2([mty',mtx'],[double(yutm),double(xutm)],'euclidean','smallest',1);
% find min and max of actual matched mop numbers
ApproxMop=mopnum(NearIdx);
%fprintf('Mop Match Range: %i to %i\n',min(ApproxMop),max(ApproxMop))
%  Get max distance of the data points from their nearest back beach point
%  Use this when defining Mop transect lines as series of closely spaced
%  points
MaxDist=max(ceil(dp)); 
%fprintf('Max Distance from Back Beach Points: %i\n',MaxDist)
%-----------------------------------------------------------------
% Now loop through approximately matched Mop numbers and do a more
% exact matching of data points to the Mop transect lines

%fprintf('Making more precise data match to each Mop transect line...\n')

N=ApproxMop;  % initialize N as the approximate Mop numbers

% Step up the coast for each Mop number. Get the survey points with appox Mop
%  numbers equal to this Mop and the next upcoast Mop and find which data
%  points are actually closest to the Mop's transect line rather than just
%  the back beach point.

for mn=min(ApproxMop)-1:max(ApproxMop)+1
    
    % mn is the target Mop number for exact nearest data matching
    if mn == min(ApproxMop)-1
%     ll=fprintf('Finding points closest to Mop %i of %i ',mn,max(ApproxMop)+1);
    else
%     fprintf(repmat('\b',1,ll))
%     ll=fprintf('Finding points closest to Mop %i of %i ',mn,max(ApproxMop)+1);
    end
    
    if mn > 0 && mn < 11594 % can't go below Mop #1
        
    % get indexes of all nearby data points with approx nearest mop number
    % equal to mn or the next upcoast mop.   
    idx=find(N >= mn & N <= mn+1); 
    
    pt=[xutm(idx) yutm(idx) xutm(idx)*0]; % make 3d [x y 0] nearby point array 

    % find distance of points to Mop mn transect line
    Mopnum=mn;
    v1=[Mop.BackXutm(Mopnum),Mop.BackYutm(Mopnum),0];
    v2=[Mop.OffXutm(Mopnum),Mop.OffYutm(Mopnum),0];
    d1 = point_to_line(pt, v1, v2); % distance to mn line
    
    % find distance of points to Mop mn+1 transect line
    Mopnum=mn+1;  
    v1=[Mop.BackXutm(Mopnum),Mop.BackYutm(Mopnum),0];
    v2=[Mop.OffXutm(Mopnum),Mop.OffYutm(Mopnum),0];
    d2 = point_to_line(pt, v1, v2); % distance to mn+1 line
    
    % Correct the overall nearest mop number vector N
    N(idx(d1 <= d2))=mn; % points closest to target mop
    N(idx(d1 > d2))=mn+1; % points closer to target mop+1
    end     
end

fprintf('Nearest Mop range: %i %i\n',min(N),max(N))
%fprintf('\n');
%fprintf('\n%i Data Pts Corrected from approx to actual nearest Mop transect.\n',length(find((ApproxMop - N) > 0)))
end
% 
% noff=0;ngood=0;
% for n=1:size(Survey,2)
% mdx=find(~ismember(Survey(n).NearestMops,Survey(n).MopStart:Survey(n).MopEnd));
% if ~isempty(mdx)
%     %fprintf('%s\n %s\n',Survey(n).File,num2str(Survey(n).NearestMops(mdx)))
%     noff=noff+1;
% else
%     ngood=ngood+1;
% end
% end
