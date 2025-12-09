% This version goes through SA surveys and regrids same truck
% lidar survey that may have been combined in previous grid
% file updates



% Plots gridded survey data points for the input MopID name or number.
%  - Finds the closest survey date to the input datenumber
%  - PlotMethod is '2Dcolor' or '3Dcolor'
%  - GridMethod is 'Delauney' 
%  - Uses data from 5 mop areas (the surrounding +/- 2 mop point
%  reach)

mpath='/volumes/group/MOPS/';

% load struct array list of moved survey data files if it exists
if exist('MovedSurveyFiles.mat','file')
    load MovedSurveyFiles.mat
end

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

for MopID=571:594
%
    fprintf('----------------- \n')
    fprintf('Mop %i \n',MopID)
    fprintf('----------------- \n')
    
% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% load selected mop's MopSurvey struct array

matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
if exist(matfile,'file') == 2
load(matfile,'SA');

for n=1:size(SA,2)
    idx=find(strcmpi({Survey.File},SA(n).File) == 1);
    SA(n).Source=Survey(idx).Source;
    SA(n).Bytes=Survey(idx).Bytes;
    SA(n).UTMzone='11 S';  
end
save(matfile,'SA');

% load exisiting grid SG struct array as old array OG
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
if exist(matfile,'file') == 2
 load(matfile,'SG');
 
for n=1:size(SG,2)
    idx=find(strcmpi({Survey.File},SG(n).File) == 1);
    SG(n).Source=Survey(idx).Source;
    SG(n).Bytes=Survey(idx).Bytes;
    SG(n).UTMzone='11 S';  
end
 
  cflag=0; % struct array change flag initialed to no change
 % remove any duplicate surveys based on the filename
 [~, idx] = unique({SG.File}, 'stable');
 if numel(idx) < size(SG,2)
     fprintf(' Removing %i Duplicate Survey File Entries from SG...\n',...
         size(SG,2)-numel(idx))    
      SG=SG(idx);
      cflag=1; % struct array change flag
 end


% first loop through SA file names and update any previously
%  found file name or file location changes
if exist('Moved','var')
    fprintf('Comparing SG to moved file list...\n')
    nmoved=1;
    mpass=0;
   while nmoved > 0 
    nmoved=0;
    mpass=mpass+1;
    for n=1:size(SG,2)
      idxm=find(strcmpi({Moved.OldName},SG(n).File) == 1);
      if ~isempty(idxm)
       nmoved=nmoved+1;
       SG(n).File=Moved(idxm(end)).NewName;
        cflag=1; % struct array change flag
       %SG(n).Bytes=Moved(idxm(end)).Bytes;
      end
    end
    fprintf('Moved %i filenames in pass %i.\n',nmoved,mpass)
   end % end while
end

% second, see if there are file name in the SG file that
%  are not found in the SA struct array (ie. reefbreak1
%  database file names have changed or been deleted for some reason)
idx0=find(~ismember(lower({SG.File}),lower({SA.File})));

%fprintf('-----------------------------------------------------\n')
fprintf('Removing data for %i filenames in SG that are not in SA.\n',numel(idx0))
%fprintf('-----------------------------------------------------\n')

% loop through from largest index to smallest and delete SG entries
% not in SA

for n=fliplr(idx0)
    SG(n)=[];
     cflag=1; % struct array change flag
end
fprintf('Size of SG: %i\n',size(SG,2))
fprintf('Size of SA: %i\n',size(SA,2))

% now find filenames in the SA struct array that
%  are not found in the SG struct array 
idx1=find(~ismember(lower({SA.File}),lower({SG.File})));

else
    SG=[];
    idx1=1:size(SA,2);
end

end

% for regridding purposes, check all the AS files
idx1=1:size(SA,2);

%fprintf('-----------------------------------------------------\n')
fprintf('There are %i new filenames in SA that are not in SG.\n',numel(idx1))
%fprintf('-----------------------------------------------------\n')

if numel(idx1) > 0
% Grid data in new SA surveys and add to SG
% 
% load neighboring 2 mops on either side
N1=[];N2=[];N3=[];N4=[];
matfile=[mpath 'M' num2str(MopNumber+2,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N1=load(matfile,'SA');end
matfile=[mpath 'M' num2str(MopNumber+1,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N2=load(matfile,'SA');end
matfile=[mpath 'M' num2str(MopNumber-1,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N3=load(matfile,'SA');end
matfile=[mpath 'M' num2str(MopNumber-2,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N4=load(matfile,'SA');end

% % Main loop through Mop survey dates
nsurv=0;

for SurvNum=idx1
    nsurv=nsurv+1;
    ngrid=SurvNum;
    
    %  see if multi single day survey by the truck
        tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
        strcmpi({SA.Source}, SA(SurvNum).Source)==1 & ...
        strcmpi({SA.Source}, 'Trk') == 1);
     %fprintf('%i %i\n',SurvNum,numel(tidx))
    
    if numel(tidx) > 1
    
    fprintf('Regridding %i %s\n',SurvNum,datestr(SA(SurvNum).Datenum))
    SG(ngrid).Mopnum=MopNumber;
    SG(ngrid).Datenum=SA(SurvNum).Datenum;
    SG(ngrid).Source=SA(SurvNum).Source;
    SG(ngrid).File=SA(SurvNum).File;
    SG(ngrid).UTMzone=SA(SurvNum).UTMzone;
     cflag=1; % struct array change flag
    
%     
      SurveyDatenum=SA(SurvNum).Datenum;
      fprintf('Mop % Survey %i of %i : %s\n',...
          MopID,nsurv,numel(idx1),datestr(SurveyDatenum))
%     
    % start aggregate x,y, z arrays      
	 Xr=[];Yr=[];Zr=[];Cr=[];
    % if multiple sets of points for same date and source
    % OTHER THAN truck lidar (Trk)
    % (Mops at boundaries of airborne data file reaches)
    %  combine for gridding
%         tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
%         strcmpi({SA.Source}, SA(SurvNum).Source)==1 & ...
%         strcmpi({SA.Source}, 'Trk') == 0); 
     tidx=SurvNum; % for regridding just use single survey     

       if ~isempty(tidx)         
	for ntidx=tidx
         Xr=[Xr' SA(ntidx).X']';
         Yr=[Yr' SA(ntidx).Y']';
         Zr=[Zr' SA(ntidx).Z']';
         Cr=[Cr' SA(ntidx).Class']';
	end
       end
      
    % inner loop through nearby Mop surveys
    
      if ~isempty(N1)
       % if it has a matching survey date, include its x,y,z points
       tidx=find(vertcat(N1.SA.Datenum) == SurveyDatenum);
       if ~isempty(tidx)         
	    for ntidx=tidx'
              Xr=[Xr' N1.SA(ntidx).X']';
              Yr=[Yr' N1.SA(ntidx).Y']';
              Zr=[Zr' N1.SA(ntidx).Z']';
              Cr=[Cr' N1.SA(ntidx).Class']';
	    end
       end
      end
      
      if ~isempty(N2)
       tidx=find(vertcat(N2.SA.Datenum) == SurveyDatenum);
       if ~isempty(tidx)         
	    for ntidx=tidx'
              Xr=[Xr' N2.SA(ntidx).X']';
              Yr=[Yr' N2.SA(ntidx).Y']';
              Zr=[Zr' N2.SA(ntidx).Z']';
              Cr=[Cr' N2.SA(ntidx).Class']';
	    end
       end
      end
      
      if ~isempty(N3)
       tidx=find(vertcat(N3.SA.Datenum) == SurveyDatenum);
       if ~isempty(tidx)         
	    for ntidx=tidx'
              Xr=[Xr' N3.SA(ntidx).X']';
              Yr=[Yr' N3.SA(ntidx).Y']';
              Zr=[Zr' N3.SA(ntidx).Z']'; 
              Cr=[Cr' N3.SA(ntidx).Class']';
	    end
       end 
      end
      
      if ~isempty(N4)
       tidx=find(vertcat(N4.SA.Datenum) == SurveyDatenum);
       if ~isempty(tidx)         
	    for ntidx=tidx'
              Xr=[Xr' N4.SA(ntidx).X']';
              Yr=[Yr' N4.SA(ntidx).Y']';
              Zr=[Zr' N4.SA(ntidx).Z']';
              Cr=[Cr' N4.SA(ntidx).Class']';
	    end
       end 
      end
      
    % end nearby survey inner loop
    
% grid the aggregated data
    
     fprintf('Gridding %i points...\n',length(Xr))
     XGutm=[];YGutm=[];Vzg=[];
     if numel(Xr) > 0
%      min(Yr)
%      max(Yr)
% bound the data gaps with NaNs to avoid gridding them
     MaxGap=250;
     %Xr(Yr > 3650000)=[];Zr(Yr > 3650000)=[];Yr(Yr > 3650000)=[];
     [x,y,z]=addNoDataAreaPoints(Xr,Yr,Zr,MaxGap);     
    
% Grid the elevation survey using Delaunay tesselation 
    zg=griddata(double(x),double(y),double(z),...
        double(min(x):max(x)),double(min(y):max(y))');
    
    [YG,XG]=find(~isnan(zg)); % find valid grid data x,y vectors
     Vzg=zg(~isnan(zg(:))); % valid zg data vector
     XGutm=min(x)-1+XG; % adjust back to utm coords
     YGutm=min(y)-1+YG;
     
     % find nearest Mop transects to each valid gridded point
     Nmop=FindNearestMopTransectsUTM(XGutm,YGutm);
     
     % reduce to valid gridded points nearest to MopNumber
     XGutm=XGutm(Nmop == MopNumber);
     YGutm=YGutm(Nmop == MopNumber);
     Vzg=Vzg(Nmop == MopNumber);
     
     end
     % add to the MopGrid struct array
      
     SG(ngrid).X=XGutm;
     SG(ngrid).Y=YGutm;
     SG(ngrid).Z=Vzg;
     
     % intialize grid point classification as unknown
     SG(ngrid).Class=zeros(size(XGutm));
     
     % replace grid point class with survey data point class at shared
     %   x,y locations.  No attempt is made to classify grid points 
     %   that were not directly surveyed.
     
     if numel(XGutm) > 0
     % use ismember to find shared 1m survey / 1m grid x,y locations
     [l,li]=ismember([Xr Yr],[XGutm YGutm],'rows');
     SG(ngrid).Class(li(l > 0))=Cr(l(l > 0));
     end
    end 
end
 
end

if cflag == 1
    
fprintf('Writing changed SG sruct array to matfile...\n')

% sort by date 
T=struct2table(SG); 
sortedT = sortrows(T, 'Datenum');
SG=table2struct(sortedT)';

% Save MopGrid struct array in its own mat file
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
save(matfile,'SG');

end

end

% % find dates with more than 1 survey file
% [~, idx] = unique([SG.Datenum], 'stable');
% DupDates=find(~ismember([1:size(SA,2)],idx));
% [~, idx2] = unique(DupDates, 'stable');
% 
% for n=unique(DupDates)
%     fprintf('%s \n',datestr(SG(n).Datenum)) 
%     idp=find([SG.Datenum] == SG(n).Datenum);
%     for i=idp
%       fprintf('%i  : %s\n',numel(SG(i).X),SA(i).File)  
%     end
% end

% % list number of SG grid points in same source duplicate date entries
% for SurvNum=1:size(SA,2)
% %   tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
% %         strcmpi({SA.Source}, SA(SurvNum).Source) == 1 &  strcmpi({SA.Source}, 'Trk') == 1); 
%   % option to include truck surveys
%   tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
%        strcmpi({SA.Source}, SA(SurvNum).Source) == 1); 
%     if numel(tidx) > 1
%      fprintf('%s \n',datestr(SA(SurvNum).Datenum)) 
%      for i=tidx
%       fprintf(' %i : %s\n',numel(SA(i).X),SA(i).File)  
%      end
%     end
% end
%     
