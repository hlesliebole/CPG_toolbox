% Adds a Single New Survey to the corresponding Mop SG matfiles
%

fprintf('UpdateSingleSurveySGmatfiles.m \n')

% SingleDate (a survey datenum) and SingleSource (survey source name)
% must be defined prio to runnnig the script.  Designed to be
% called by UpdateSingleSurveySASGSM.m

load SurveyMasterList.mat

idx=find([Survey.Datenum] == SingleDate & contains({Survey.Source},SingleSource) == 1);

load('MopTableUTM.mat','Mop');  % Load "Mop" table array


if numel(idx) == 1
    
Mop1=Survey(idx).MopStart;
Mop2=Survey(idx).MopEnd;
%Mop1=585;Mop2=585;

for MopID=Mop1:Mop2

    fprintf('Mop %i \n',MopID)
    
% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% load selected mop's MopSurvey struct array

matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
if exist(matfile,'file') == 2
load(matfile,'SA');

SAidx=find([SA.Datenum] == datenum(2021,10,7) & contains({SA.Source},'Gps') == 1);
if numel(SAidx) == 1
    
% load exisiting grid SG struct array as old array OG
matfile=['M' num2str(MopNumber,'%5.5i') 'SG.mat'];
if exist(matfile,'file') == 2
 load(matfile,'SG');
 OG=SG;
else
 clear OG
 OG(1).Datenum=0;
 OG(2).Source=' ';
end

% SG array should be one survey smaller than SA. If
%  not, stop.

if size(SA,2) - size(SG,2) == 1

% Initializing new SG struct array up to the
%  survey idx-1 to match SA array
SG=OG(1:SAidx-1);

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

% New survey date
SurvNum=SAidx;
  
      SurveyDatenum=SA(SurvNum).Datenum;
      fprintf('Adding Survey %i of %i : %s\n',...
          SurvNum,size(SA,2),datestr(SurveyDatenum))
    
    % start aggregate x,y, z arrays      
	 Xr=[];Yr=[];Zr=[];Cr=[];
    % if multiple sets of points for same date and source 
    % (Mops at boundaries of airborne data file reaches)
    %  combine for gridding
        tidx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
        strcmpi({SA.Source}, SA(SurvNum).Source)==1); 
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
     
     SG(SurvNum).Mopnum=SA(SurvNum).Mopnum;
     SG(SurvNum).Datenum=SA(SurvNum).Datenum;
     SG(SurvNum).Source=SA(SurvNum).Source;
     SG(SurvNum).File=SA(SurvNum).File;
     SG(SurvNum).UTMzone=SA(SurvNum).UTMzone;
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
 

% add any existing surveys after the added one
if size(OG,2)+1 > size(SG,2)
    SG(SAidx+1:size(OG,2)+1)=OG(SAidx:end);
end

% Save MopGrid struct array in its own mat file
matfile=['M' num2str(MopNumber,'%5.5i') 'SG.mat'];
save(matfile,'SG');

else
    fprintf('SA is not 1 survey larger than SG. Skipping.\n')
end

else
   fprintf('SA does not contain the survey. Skipping.\n') 
end

else
   fprintf('SA mat files does not exist. Skipping.\n') 
end 

end  % end mop loop

else
    fprintf('No Matching Survey for that date and source.\n')
end