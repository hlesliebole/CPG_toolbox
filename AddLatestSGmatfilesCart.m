%function SG=MakeMopGrid(MopID)


% Plots gridded survey data points for the input MopID name or number.
%  - Finds the closest survey date to the input datenumber
%  - PlotMethod is '2Dcolor' or '3Dcolor'
%  - GridMethod is 'Delauney' 
%  - Uses data from 5 mop areas (the surrounding +/- 2 mop point reach)
% folder containing CPG Mop mat files

MopPath='/volumes/group/MOPS/';

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

for MopID=578:584

    fprintf('%i \n',MopID)
    
% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% load selected mop's MopSurvey struct array

matfile=[MopPath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
%if exist(matfile,'file') == 2
load(matfile,'SA');

% load exisiting grid SG struct array as old array OG
matfile=[MopPath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
load(matfile,'SG');

% only add latest "ns" SA surveys if SG has ns less than SA
%ns=1;
%if size(SA,2)-size(SG,2) == ns
ns=size(SA,2)-size(SG,2);
if ns < 42

% load neighboring 2 mops on either side
N1=[];N2=[];N3=[];N4=[];
matfile=['M' num2str(MopNumber+2,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N1=load(matfile,'SA');end
matfile=['M' num2str(MopNumber+1,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N2=load(matfile,'SA');end
matfile=['M' num2str(MopNumber-1,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N3=load(matfile,'SA');end
matfile=['M' num2str(MopNumber-2,'%5.5i') 'SA.mat'];
if exist(matfile,'file');N4=load(matfile,'SA');end

% Process last ns SA survey dates
for SurvNum=size(SA,2)-ns+1:size(SA,2)
%for SurvNum=1:1
    
      SurveyDatenum=SA(SurvNum).Datenum;
      fprintf('Survey %i of %i : %s\n',...
          SurvNum,size(SA,2),datestr(SurveyDatenum))
   
    % Survey date and source already exist in old grid struct
    %  area, reassign to new SG struct array
%     old=find( [OG.Datenum] == SA(SurvNum).Datenum & ...
%         strcmpi({OG.Source}, SA(SurvNum).Source)==1);
%     if isempty(old)  % new survey data, grid it and add to SG struct array
%         fprintf(' New Data %s %s\n',datestr(SA(SurvNum).Datenum),SA(SurvNum).Source)
%         if numel(SA(SurvNum).X) >= 3 
    
    % start aggregate x,y, z arrays      
	 Xr=[];Yr=[];Zr=[];Cr=[];
%     % if multiple sets of points for same date and source 
%     % (Mops at boundaries of airborne data file reaches)
%     %  combine for gridding
%         tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
%         strcmpi({SA.Source}, SA(SurvNum).Source)==1); 
%        if ~isempty(tidx)         
% 	for ntidx=tidx
         ntidx=SurvNum; % don't look for multiples
         Xr=[Xr' SA(ntidx).X']';
         Yr=[Yr' SA(ntidx).Y']';
         Zr=[Zr' SA(ntidx).Z']';
         Cr=[Cr' SA(ntidx).Class']';
%	end
%       end
      
    % inner loop through nearby Mop surveys
    
      if ~isempty(N1)
       % if it has a matching survey date, include its x,y,z points
       %tidx=find(vertcat(N1.SA.Datenum) == SurveyDatenum);
       tidx=find(strcmpi({N1.SA.File}, SA(SurvNum).File));
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
       %tidx=find(vertcat(N2.SA.Datenum) == SurveyDatenum);
       tidx=find(strcmpi({N2.SA.File}, SA(SurvNum).File));
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
       %tidx=find(vertcat(N3.SA.Datenum) == SurveyDatenum);
       tidx=find(strcmpi({N3.SA.File}, SA(SurvNum).File));
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
       %tidx=find(vertcat(N4.SA.Datenum) == SurveyDatenum);
       tidx=find(strcmpi({N4.SA.File}, SA(SurvNum).File));
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
     
     % add to the MopGrid struct array
     
     % copy meta info from SA struct array
     SG(SurvNum).Mopnum=SA(SurvNum).Mopnum;
     SG(SurvNum).Datenum=SA(SurvNum).Datenum;
     SG(SurvNum).Source=SA(SurvNum).Source;
     SG(SurvNum).File=SA(SurvNum).File;
     SG(SurvNum).UTMzone=SA(SurvNum).UTMzone;
     
     % gridded data points
     SG(SurvNum).X=XGutm;
     SG(SurvNum).Y=YGutm;
     SG(SurvNum).Z=Vzg;
     
     % intialize grid point classification as unknown
     SG(SurvNum).Class=zeros(size(XGutm));
     
     % replace grid point class with survey data point class at shared
     %   x,y locations.  No attempt is made to classify grid points 
     %   that were not directly surveyed.
     
     % use ismember to find shared 1m survey / 1m grid x,y locations
     [l,li]=ismember([Xr Yr],[XGutm YGutm],'rows');
     SG(SurvNum).Class(li(l > 0))=Cr(l(l > 0));
 
%end

% Save MopGrid struct array in its own mat file
matfile=[MopPath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
save(matfile,'SG');

end

else
    fprintf('SA size %i not 1 larger than SG size %i\n',size(SA,2),size(SG,2))
end

end 
