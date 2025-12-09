clearvars
CpgDefineMopPath
%mpath='/volumes/group/mops/';

load SurveyMasterListWithMops.mat
load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% find all the Mop numbers with survey data
MopNums=unique([Survey.NearestMops]);
MopNums=sort(MopNums);

for MopID=MopNums
    if MopID < 1690 | (MopID > 2650 & MopID < 2701)
    m=MopID;
    %fprintf('%i\n',m)
ndx=find(cellfun(@(v)any(v(:) == m),{Survey.NearestMops}));
% fprintf('The current Survey struct array has %i surveys for this Mop.\n',...
%     numel(ndx))
    
%--------

% first just read in the smaller SMfiles inventory struct for this mop 
%  to check for changes

%load([ mpath 'M'  num2str( m , '%5.5i' )  'SA.mat' ],'SAfiles');
% if size(SAfiles,2) ~= numel(ndx)
%     fprintf('Mop %i SA and Master Survey mismatch: %i %i\n',...
%         m,size(SAfiles,2),numel(ndx))
% end
 
% read in the SM file if it exists

SGfiles=[];
if exist([ mpath 'M'  num2str( m , '%5.5i' )  'SG.mat' ],'file')
load([ mpath 'M'  num2str( m , '%5.5i' )  'SG.mat' ],'SGfiles');
end
% if size(SGfiles,2) ~= numel(ndx)
%     fprintf('Mop %i SG and Master Survey mismatch: %i %i\n',...
%         m,size(SGfiles,2),numel(ndx))
% end

   if size(SGfiles,2) > 0
   % combine filedatenums and filenames into single strings for comparison
   % of both the file names and creation dates at the same time
   Gdf=strcat(string([SGfiles.FileDatenum]),{SGfiles.File}); % SG string
   Sdf=strcat(string([Survey(ndx).FileDatenum]),{Survey(ndx).File}); % Survey string
   %kdx=find(~ismember({SGfiles.File},{Survey(ndx).File}))
   ndel=find(~ismember(Gdf,Sdf)); % SG file that is not a member of Survey
   nadd=find(~ismember(Sdf,Gdf)); % Survey file not a member of SG
   if numel(ndel) > 0
       fprintf(' Mop %i Deleting %i surveys\n',m,numel(ndel))
   end
   if numel(nadd) > 0
       fprintf(' Mop %i Adding %i surveys\n',m,numel(nadd))
   end
   
   else
       fprintf('No SG data for Mop %i.  %i surveys being added.\n',m,numel(ndx))
   end
   
if numel(ndel) > 0 || numel(nadd) > 0
 load([ mpath 'M'  num2str( m , '%5.5i' )  'SG.mat' ],'SG');
end

if exist('SG','var')
if ~isempty(SG)
if numel(ndel) > 0 & ndel <= size(SG,2)
    SG(ndel)=[];
    SGfiles=struct('File',{SG.File},'FileDatenum',...
            num2cell([SG.FileDatenum]));
    save([ mpath 'M'  num2str( m , '%5.5i' )  'SG.mat' ],'SG','SGfiles');
end
end
end
   
if numel(nadd) > 0
    
MopNumber=m;

% load SA array for this mop
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
if exist(matfile,'file');load(matfile,'SA');end

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

% Main loop through new Mop survey dates
nd=size(SG,2);

   for SurveyNum=ndx(nadd)
      fprintf('Adding: %s\n',Survey(SurveyNum).File)
      % find SA survey number matching new Survey file name
        
%       SurveyDatenum=SA(SurvNum).Datenum;
%       fprintf('Mop %i Survey %i of %i : %s \n%s\n',...
%           MopID,SurvNum,size(SA,2),datestr(SurveyDatenum),SA(SurvNum).File)
%       
    % start aggregate x,y, z arrays 
         Xr=[];Yr=[];Zr=[];Cr=[];
        
         SurvNum=find(strcmpi({SA.File},Survey(SurveyNum).File));
        if ~isempty(SurvNum) 
         [x,y,z,c]=CombineSurveyDateSourceData(SA,SurvNum);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
        end
      
    % inner loop through nearby Mop surveys
    %
    %  Only include N1 and N4 (+/- 2 mops away) if it is a "Gps" survey
    
%     if strcmpi(Survey(SurveyNum).Source,'Gps') == 1     
      if ~isempty(N1)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N1.SA.File},Survey(SurveyNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N1.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
      end
%     end
      
       if ~isempty(N2)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N2.SA.File},Survey(SurveyNum).File) == 1);
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
       tidx=find(strcmp({N3.SA.File},Survey(SurveyNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N3.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
      end
       
%     if strcmpi(Survey(SurveyNum).Source,'Gps') == 1    
       if ~isempty(N4)
       % if it has a matching survey file, include its x,y,z points
       tidx=find(strcmp({N4.SA.File},Survey(SurveyNum).File) == 1);
       if ~isempty(tidx) 
             [x,y,z,c]=CombineSurveyDateSourceData(N4.SA,tidx);
              Xr=[Xr' x']';
              Yr=[Yr' y']';
              Zr=[Zr' z']';
              Cr=[Cr' c']';
       end
       end 
%     end
      
    % end nearby survey inner loop
    
    if ~isempty(Xr)
        
    % remove any points flagged as bad with a .Class=-999
    Xr(Cr == -999)=[];Yr(Cr == -999)=[];Zr(Cr == -999)=[];
    
% % first grid the aggregated data in alongshore  Mop coordinate
% %    frame
%      fprintf('Gridding %i points in Mop alongshore coords...\n',numel(Xr))
%      [Xa,Ya,Za]=AlongshoreGrid(Xr,Yr,Zr);
  
     % option to skip initial alongshore gridding and just do normal
     %  gridding
     Xa=Xr;Ya=Yr;Za=Zr;

     % check for duplicate x, y, z
     [uxy,ia,ic]=unique([Xa,Ya],'rows');uniqxyz=size(uxy,1);
     % remove those points
     Xa=Xa(ia);Ya=Ya(ia);Za=Za(ia);
    
     fprintf('Gridding %i points in UTM coords...\n',numel(Xa))
     
     
% bound the data gaps with NaNs to avoid gridding them
     MaxGap=250;
     %Xr(Yr > 3650000)=[];Zr(Yr > 3650000)=[];Yr(Yr > 3650000)=[];
     if numel(Xr) > 5  % only grid if there are sufficient points
     [x,y,z]=addNoDataAreaPoints(Xa,Ya,Za,MaxGap);  
     
     % check for duplicate x, y, z
     [uxy,ia,ic]=unique([x,y],'rows');
     % remove those points
     x=x(ia);y=y(ia);z=z(ia);   
    
% Grid the elevation survey using Delaunay tesselation 
    zg=griddata(double(x),double(y),double(z),...
        double(min(x):max(x)),double(min(y):max(y))');
    
    [YG,XG]=find(~isnan(zg)); % find valid grid data x,y vectors
     Vzg=zg(~isnan(zg(:))); % valid zg data vector
     XGutm=min(x)-1+XG; % adjust back to utm coords
     YGutm=min(y)-1+YG;
     
     if ~isempty(XGutm)
     
     % find nearest Mop transects to each valid gridded point
     %Nmop=FindNearestMopTransectsUTM(XGutm,YGutm);
     Nmop=XY2MopNumsV2(XGutm,YGutm,Mop);
     
     % reduce to valid gridded points nearest to MopNumber
     XGutm=XGutm(Nmop == MopNumber);
     YGutm=YGutm(Nmop == MopNumber);
     Vzg=Vzg(Nmop == MopNumber);
     
     % add to the MopGrid struct array
     nd=nd+1;
     SG(nd).Mopnum=MopNumber;
     SG(nd).UTMzone=SA(1).UTMzone;
     SG(nd).File=Survey(SurveyNum).File;
     SG(nd).FileDatenum=Survey(SurveyNum).FileDatenum;
     SG(nd).Source=Survey(SurveyNum).Source;
     SG(nd).Datenum=Survey(SurveyNum).Datenum; 
     SG(nd).X=XGutm;
     SG(nd).Y=YGutm;
     SG(nd).Z=Vzg;
     
     % intialize grid point classification as unknown
     SG(nd).Class=zeros(size(XGutm));
     
     % replace grid point class with survey data point class at shared
     %   x,y locations.  No attempt is made to classify grid points 
     %   that were not directly surveyed.
     
     % use ismember to find shared 1m survey / 1m grid x,y locations
     [l,li]=ismember([Xr Yr],[XGutm YGutm],'rows');
     SG(nd).Class(li(l > 0))=Cr(l(l > 0));
     
     else
      fprintf('No data for this Mop\n'); 
     end
     
     else
      fprintf('No data for this Mop\n');
     end
     
   else
      fprintf('No SA data for this Mop\n');       
   end



   end
   
% % sort survey list by survey date
    if size(SG,2) > 1
     T=struct2table(SG);
     sortedT = sortrows(T, 'Datenum');
     SG=table2struct(sortedT)';
    
    
    % remove any duplicate surveys based on the filename
     [~, idx] = unique({SG.File}, 'stable');
     if numel(idx) < size(SG,2)
         fprintf(' Removing %i Duplicate Survey File Entries from SA...\n',...
             size(SG,2)-numel(idx)) 
          ddx=find(~ismember(1:size(SG,2),idx));
          for n=ddx
              fprintf('%s\n',SG(n).File)
          end
          SG=SG(idx);
     else
         %fprintf('No duplicate file entries found.\n')
     end
        
    end
   
% Save MopGrid struct array in its own mat file
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat'];
if size(SG,2) > 0
SGfiles=struct('File',{SG.File},'FileDatenum',...
            num2cell([SG.FileDatenum]));
else
    SG=[];
    SGfiles=[];
end

save(matfile,'SG','SGfiles');

end
    end
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


