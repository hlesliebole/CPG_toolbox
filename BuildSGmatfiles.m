%function SG=MakeMopGrid(MopID)

%mpath='/volumes/group/MOPS/';
clearvars
CpgDefineMopPath

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

Mop1=560;
Mop2=600;
for MopID=Mop1:Mop2
    fprintf('%i of %i to %i \n',MopID,Mop1,Mop2)
    
% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% load selected mop's MopSurvey struct array

matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');

% Initializing new SG struct array using the SA
%  struct array. The x, z, data fields will be replaced with 
%  1m gridded data points.
SG=SA;

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

% Main loop through Mop survey dates
for SurvNum=1:size(SA,2)
%for SurvNum=1:1
      SurveyDatenum=SA(SurvNum).Datenum;
      fprintf('Survey %i of %i : %s \n%s\n',...
          SurvNum,size(SA,2),datestr(SurveyDatenum),SA(SurvNum).File)
      
    % start aggregate x,y, z arrays 
         [Xr,Yr,Zr,Cr]=CombineSurveyDateSourceData(SA,SurvNum);
      
    % inner loop through nearby Mop surveys
    %
    %  Only include N1 and N4 (+/- 2 mops away) if it is a "Gps" survey
    
    if strcmpi(SA(SurvNum).Source,'Gps') == 1      
      if ~isempty(N1)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N1.SA.File},SA(SurvNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N1.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
      end
    end
      
       if ~isempty(N2)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N2.SA.File},SA(SurvNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N2.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
      end
      
       if ~isempty(N3)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N3.SA.File},SA(SurvNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N3.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
      end
       
    if strcmpi(SA(SurvNum).Source,'Gps') == 1    
       if ~isempty(N4)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N4.SA.File},SA(SurvNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N4.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
       end 
    end
      
    % end nearby survey inner loop
    
% grid the aggregated data
    
     fprintf('Gridding %i points...\n',numel(Xr))
     
     if numel(Xr) > 5  % only grid if there are sufficient points
% bound the data gaps with NaNs to avoid gridding them
     MaxGap=250;
     Xr(Yr > 3650000)=[];Zr(Yr > 3650000)=[];Yr(Yr > 3650000)=[];
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
     
     else
      SG(SurvNum).X=[];
      SG(SurvNum).Y=[];
      SG(SurvNum).Z=[];
      SG(SurvNum).Class=[];
     end
     
end

% Save MopGrid struct array in its own mat file
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
save(matfile,'SG');

end

 