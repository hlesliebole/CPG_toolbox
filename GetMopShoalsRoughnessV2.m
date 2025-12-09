function RS=GetMopShoalsRoughnessV2(MopNumber,WinSize,MinWinPts)

%% eg. RS=GetMopShoalsRoughnessV2(667,11,5) % for cardiff seaside reef

%% input:
%
% MopNumber = single Mop number to process
% WinSize = 2d x,y spatial filter window size in meters. Odd numbered size is better.
% MinWinPts = Minimum number of points within the moving 2d window to return a
%             a roughness (sigma) value for the center location.

%%  returns "roughness" grid struct array RS of SHOALS spatially-filtered standard deviations 
%      RS(N).SurveyDatenum
%      RS(N).SpatialWindowSize;
%      RS(N).SpatialMinPoints;
%      RS(N).Xutm 
%      RS(N).Yutm  
%      RS(N).Z % spatial filter avg elev NAVD88
%      RS(N).Sigma   % spatial filter standard deviation "roughness"
%
%  N = 1 = Shoals 2009 survey
%  N = 2 = SHoals 2014 survey
%  N = 3 = Combined 2009 and 2014 surveys, RS(3).SurveyDatenum set to 0.

%% dependencies
%
%  uses MOPS/toolbox m-scripts:
%
%     UTM2MopxshoreX.m
%       - FindNearestMopTransectsUTM.m
%       - GetTransectLines.m
%
%     XY2MopNumsV2.m
%
%  and the MOPS/toolbox .mat file
%
%     MopTableUTM.mat
%
load('MopTableUTM.mat','Mop');

for n=1:3
RS(n).SpatialWindowSize=WinSize;
RS(n).SpatialMinPoints=MinWinPts;
RS(n).SurveyDatenum=[];
RS(n).Xutm=[];
RS(n).Yutm=[];
RS(n).Z=[];
RS(n).Sigma=[];
end

%% Combine 1m avg Mop SA struct arrays +/- either side of Mop 
CS=SAcombineMops(MopNumber-1,MopNumber+1);

%% get global x,y limits for the composite grid
%   use all of the shoals surveys in both years 
sidx=find(strcmp({CS.Source},'USACE'));
xmin=[];xmax=[];ymin=[];ymax=[];
for n=sidx  % loop through shoals surveys
    
  % % reduce to points with positive MOP shoreline xshore locations
  % %   to remove any estuary points
 [Nmop,xs]=UTM2MopxshoreX([CS(n).X],[CS(n).Y]);
  ldx=find(xs < -50);
  %fprintf('%i %i\n',n,numel(ldx))
  CS(n).X(ldx)=[];
  CS(n).Y(ldx)=[];
  CS(n).Z(ldx)=[];
  
  % add to min max calcs
 xmin=min([xmin [CS(n).X]']);xmax=max([xmax [CS(n).X]']);
 ymin=min([ymin [CS(n).Y]']);ymax=max([ymax [CS(n).Y]']);
 
end

%% Process the 2009 and 2014 surveys separately and make grid points for both
%
% make 2d X,Y grid arrays
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);

%------------------
% 2009
iyr=2009;
sidx=find(strcmp({CS.Source},'USACE') & year(datetime([CS.Datenum],'convertfrom','datenum')) == iyr);%%  
%progressively overlay Shoals 1m average data onto the grid 
%  so any grid points with redundant info keep the most recent survey
Z=nan(size(X));% initialize composite elev grid as no data NaNs  
for n=sidx % loop through shoal survey dates  
    gdx=sub2ind(size(X),round(CS(n).Y)-ymin+1,round(CS(n).X)-xmin+1);% 1d grid indices with data 
    Z(gdx)=CS(n).Z; % overlay on grid
end
% get 2d spatially filtered avg elev and standard deviation (roughness)
[avgZ1,stdZ1]=runfilt2d(WinSize,Z,MinWinPts); %use WinSize and monNP settings

R1=nan(size(X));% initialize 2009 roughness grid as no data NaNs
D1=nan(size(X));% initialize 2009 depth grid as no data NaNs
idx=find(~isnan(stdZ1(:))); 
% assign them to nearest mops
if numel(idx) > 0

N=XY2MopNumsV2(X(idx),Y(idx),Mop);

if ~isempty(idx);R1(idx)=stdZ1(idx);D1(idx)=avgZ1(idx);end

RS(1).SurveyDatenum=CS(sidx(end)).Datenum;
RS(1).Xutm=X(idx(N == MopNumber));
RS(1).Yutm=Y(idx(N == MopNumber));
RS(1).Z=D1(idx(N == MopNumber));
RS(1).Sigma=R1(idx(N == MopNumber));

end

%------------------
% 2014
iyr=2014;
sidx=find(strcmp({CS.Source},'USACE') & year(datetime([CS.Datenum],'convertfrom','datenum')) == iyr);%%  
%progressively overlay Shoals 1m average data onto the grid 
%  so any grid points with redundant info keep the most recent survey
Z=nan(size(X));% initialize composite elev grid as no data NaNs  
for n=sidx % loop through shoal survey dates  
    gdx=sub2ind(size(X),round(CS(n).Y)-ymin+1,round(CS(n).X)-xmin+1);% 1d grid indices with data 
    Z(gdx)=CS(n).Z; % overlay on grid
end
% get 2d spatially filtered avg elev and standard deviation (roughness)
[avgZ2,stdZ2]=runfilt2d(WinSize,Z,MinWinPts); % use WinSize and monNP settings

R2=nan(size(X));% initialize 2014 roughness grid as no data NaNs
D2=nan(size(X));% initialize 2014 depth grid as no data NaNs
idx=find(~isnan(stdZ2(:))); 
% assign them to nearest mops
if numel(idx) > 0
N=XY2MopNumsV2(X(idx),Y(idx),Mop);
if ~isempty(idx);R2(idx)=stdZ2(idx);D2(idx)=avgZ2(idx);end

RS(2).SurveyDatenum=CS(sidx(end)).Datenum;
RS(2).Xutm=X(idx(N == MopNumber));
RS(2).Yutm=Y(idx(N == MopNumber));
RS(2).Z=D2(idx(N == MopNumber));
RS(2).Sigma=R2(idx(N == MopNumber));
end


%% blend grid data together.
%
R=nan(size(X));% initialize composite roughness grid as no data NaNs
D=nan(size(X));% initialize composite depth grid as no data NaNs
% 2009 has data, 2014 does not
idx=find(isnan(stdZ2(:)) & ~isnan(stdZ1(:))); 
% assign them to nearest mops
%N=XY2MopNumsV2(X(idx),Y(idx),Mop);
if ~isempty(idx);R(idx)=stdZ1(idx);D(idx)=avgZ1(idx);end
% 2014 has data, 2009 does not
idx=find(isnan(stdZ1(:)) & ~isnan(stdZ2(:)));
if ~isempty(idx);R(idx)=stdZ2(idx);D(idx)=avgZ2(idx);end
% both have data, 2009 has lower elev
idx=find(avgZ1(:) <= avgZ2(:)); 
if ~isempty(idx);R(idx)=stdZ1(idx);D(idx)=avgZ1(idx);end
% both have data, 2014 has lower elev
idx=find(avgZ2(:) <= avgZ1(:)); 
if ~isempty(idx);R(idx)=stdZ2(idx);D(idx)=avgZ2(idx);end

%% Just save composite roughness grid points with data
% 1d indices of grid points with data
idx=find(~isnan(R(:))); 
% assign them to nearest mops
if numel(idx) > 0
N=XY2MopNumsV2(X(idx),Y(idx),Mop);
% retain only the x,y,z,roughness data in the input mop area
RS(3).SurveyDatenum=0;
RS(3).Xutm=X(idx(N == MopNumber));
RS(3).Yutm=Y(idx(N == MopNumber));
RS(3).Z=D(idx(N == MopNumber));
RS(3).Sigma=R(idx(N == MopNumber));
end

end

