%  ReProcessSingleSurvey.m
mpath='/Volumes/group/MOPS/';

fn='/Volumes/group/topobathy/20190717_00567_00635_torreydelmar_jumbo/filtered_clean20190717.llnezts.navd88';
load SurveyMasterListWithMops.mat
idx=find(strcmpi({Survey.File},fn) == 1);

if ~isempty(idx)
    
Survey=Survey(idx); % reduce survey list to the single survey

MopStart=1; % MX border 
MopEnd=11594; % OR border

% First, cycle through all existing Mop SA files any struct array
%  entries for file name.

for m=MopStart:MopEnd
    matfile=[mpath 'M' num2str(m,'%5.5i') 'SA.mat'];
    if exist(matfile,'file')
      ll=fprintf('%s\n',matfile);
      load(matfile,'SA');
      idx=find(strcmpi({SA.File},fn) == 1);
      if ~isempty(idx)
          SA(idx)=[];
          fprintf('Removed Entry. Saving change.\n');
          save(matfile,'SA');
      else
          fprintf(repmat('\b',1,ll))
      end
    end
end

%----------------------
% now process the single file
%----------------------

nsurv(MopStart:MopEnd)=0;
% 

for n=1:1
    %fprintf('Survey %i of %i\n',n,size(Survey,2))
    
    [x,y,z,c,utmzone]=readSurveyFileUTM2(Survey(n).Source,Survey(n).File);
  if ~isempty(x)
    ext=Survey(n).File(end-2:end); % get filename extension
    if(strcmpi(ext,'tif') == 1 || strcmpi(Survey(n).Source,'UTAir') == 1)
    else
        % spatially average to 1m if *not* a tif file or UT air lidar
         [x,y,z,c]=SpatialAverageUTM(x,y,z,c,1);
    end
    [nmop]=XY2MopNumsV2(x,y,Mop);  % nearest mop area for each point
    % save in SA struct arrays for each mop
    if ~isempty(nmop)
      for m=unique(nmop)
        if m >= MopStart && m <= MopEnd 
            nsurv(m)=nsurv(m)+1;
            N(m).SA(nsurv(m)).Mopnum=m;
            N(m).SA(nsurv(m)).Datenum=Survey(n).Datenum;
            N(m).SA(nsurv(m)).Source=Survey(n).Source;
            N(m).SA(nsurv(m)).File=Survey(n).File;
            N(m).SA(nsurv(m)).UTMzone=utmzone(1,:);
            N(m).SA(nsurv(m)).X=x(nmop == m);
            N(m).SA(nsurv(m)).Y=y(nmop == m);
            N(m).SA(nsurv(m)).Z=z(nmop == m);
            N(m).SA(nsurv(m)).Class=c(nmop == m);
            N(m).SA(nsurv(m)).Bytes=Survey(n).Bytes;
        end
      end
    end
  end
end

% save SA data as individual mat files
fprintf('Saving data as SA matfiles...\n')
nsa=0; % matfile counter
for m=1:size(N,2)
  if ~isempty(N(m).SA)
    fprintf('%i\n',m)
   if exist([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'file')
    load([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'SA');
   else
    SA=[];
   end
    SA2=N(m).SA;
    SA=[SA SA2];
    % % sort survey list by survey date
    T=struct2table(SA);
    sortedT = sortrows(T, 'Datenum');
    SA=table2struct(sortedT)';
    save([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'SA');
    nsa=nsa+1;
  end
end
fprintf('%i matfiles updated.\n',nsa)

else
  fprintf('File not found in SurveyMasterListWithMops.mat Survey struct array:\n%s\n',...
      fn)  
end

%-------------------------------------------------------------------

function [xutm,yutm,zavg,cmode]=SpatialAverageUTM(xutm,yutm,z,c,Res)

%----------------------------------------
% reduce to 1m spatial averages  
%----------------------------------------
 
%Res=1; % 1m spatial resolution

% round survey x,y to desired resolution
xr=Res*round(xutm/Res); % round to Res meter spatial resolution
yr=Res*round(yutm/Res); % 

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(xr);
[uy, ~, yidx] = unique(yr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average of z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
%tavg = accumarray([xidx(:), yidx(:)], t.')./zcount;
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data

% final 1m average data vectors
xutm=ux(i);yutm=uy(j);
zavg=zavg(ii);
cmode=cmode(ii);
%tavg=tavg(ii);

%fprintf('Min-Max Classification Mode: %d %d\n',min(cmode),max(cmode))

%--------------------------------------------------------
% option to remove spatially averaged data based on less 
%  than N survey points
%
% zcount=zcount(ii); 
% N=3;
% i=find(zcount < N);xutm(i)=[];yutm(i)=[];zavg(i)=[];
%  zcount(i)=[];cmode(i)=[];
%--------------------------------------------------------

%end % end if for npts > 0

fprintf(1,...
    'Reduced to %g , %g x %g meter spatially averaged survey points\n',...
    length(xutm),Res,Res);

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

