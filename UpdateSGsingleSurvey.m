%function SG=MakeMopGrid(MopID)

MopID=581;
SurvNum=366;

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

% load exisiting grid SG struct array 
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
if exist(matfile,'file') == 2
load(matfile,'SG');

%  regrid SurvNum if file names match for this survey 
%  number in the SA and SG grids
if strcmp(SA(SurvNum).File,SG(SurvNum).File) == 1  

% Grid data in SA surveys and update SG
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

for SurvNum=366
    
    ngrid=366;
    
    SG(ngrid).Mopnum=MopNumber;
    SG(ngrid).Datenum=SA(SurvNum).Datenum;
    SG(ngrid).Source=SA(SurvNum).Source;
    SG(ngrid).File=SA(SurvNum).File;
    SG(ngrid).UTMzone=SA(SurvNum).UTMzone;
    %     
    SurveyDatenum=SA(SurvNum).Datenum;
   
%     
    % start aggregate x,y, z arrays      
	 Xr=[];Yr=[];Zr=[];Cr=[];
    
 
         Xr=[Xr' SA(SurvNum).X']';
         Yr=[Yr' SA(SurvNum).Y']';
         Zr=[Zr' SA(SurvNum).Z']';
         Cr=[Cr' SA(SurvNum).Class']';
         fprintf('%8i %s\n',numel(SA(SurvNum).Z,SA(SurvNum).File))
	
    % inner loop through nearby Mop surveys
    
      if ~isempty(N1)
      
       ntidx=find(strcmpi({N1.SA.File}, SA(SurvNum).File)==1);
       if ~isempty(ntidx)         
              Xr=[Xr' N1.SA(ntidx).X']';
              Yr=[Yr' N1.SA(ntidx).Y']';
              Zr=[Zr' N1.SA(ntidx).Z']';
              Cr=[Cr' N1.SA(ntidx).Class']';
              fprintf('%8i %s\n',numel(N1.SA(ntidx).Z,N1.SA(ntidx).File))
       else
              fprintf('%8i %s\n',0,N1.SA(ntidx).File)
       end
      end
      
      if ~isempty(N2)
       
       ntidx=find(strcmpi({N2.SA.File}, SA(SurvNum).File)==1);
       if ~isempty(ntidx)         
              Xr=[Xr' N2.SA(ntidx).X']';
              Yr=[Yr' N2.SA(ntidx).Y']';
              Zr=[Zr' N2.SA(ntidx).Z']';
              Cr=[Cr' N2.SA(ntidx).Class']';
              fprintf('%8i %s\n',numel(N2.SA(ntidx).Z,N2.SA(ntidx).File))
       else
              fprintf('%8i %s\n',0,N2.SA(ntidx).File)
       end
      end
      
      if ~isempty(N3)
      
       ntidx=find(strcmpi({N3.SA.File}, SA(SurvNum).File)==1);
       if ~isempty(ntidx)         
              Xr=[Xr' N3.SA(ntidx).X']';
              Yr=[Yr' N3.SA(ntidx).Y']';
              Zr=[Zr' N3.SA(ntidx).Z']';
              Cr=[Cr' N3.SA(ntidx).Class']';
              fprintf('%8i %s\n',numel(N3.SA(ntidx).Z,N3.SA(ntidx).File))
       else
              fprintf('%8i %s\n',0,N3.SA(ntidx).File)
       end
      end
        
      if ~isempty(N4)
      
       ntidx=find(strcmpi({N4.SA.File}, SA(SurvNum).File)==1);
       if ~isempty(ntidx)         
              Xr=[Xr' N4.SA(ntidx).X']';
              Yr=[Yr' N4.SA(ntidx).Y']';
              Zr=[Zr' N4.SA(ntidx).Z']';
              Cr=[Cr' N4.SA(ntidx).Class']';
              fprintf('%8i %s\n',numel(N4.SA(ntidx).Z,N4.SA(ntidx).File))
       else
              fprintf('%8i %s\n',0,N4.SA(ntidx).File)
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

    
fprintf('Writing changed SG sruct array to matfile...\n')

% sort by date 
T=struct2table(SG); 
sortedT = sortrows(T, 'Datenum');
SG=table2struct(sortedT)';

% Save MopGrid struct array in its own mat file
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
%save(matfile,'SG');

end

else
    fprintf('SA and SG struct arrays do not match:\n %s\n %s\n',...
        SA(SurvNum).File,SG(SurvNum).File)  
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

% list number of SG grid points in same source duplicate date entries
% for SurvNum=1:size(SA,2)
%   tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
%         strcmpi({SA.Source}, SA(SurvNum).Source) == 1 &  strcmpi({SA.Source}, 'Trk') == 0); 
%   % option to include truck surveys
%   %tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
%   %      strcmpi({SA.Source}, SA(SurvNum).Source) == 1); 
%     if numel(tidx) > 1
%     if numel(tidx) > 1
%      fprintf('%s \n',datestr(SG(SurvNum).Datenum)) 
%      for i=tidx
%       fprintf('%i  : %s\n',numel(SG(i).X),SG(i).File)  
%      end
%     end
% end
    
