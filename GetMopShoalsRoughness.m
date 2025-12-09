function RS=GetMopShoalsRoughness(MopNumber)

%%  returns "roughness" grid struct array RS of SHOALS spatially-filtered standard deviations 
%      RS.Mopnum
%      RS.SpatialWindowSize;
%      RS.SpatialMinPoints;
%      RS.Xutm 
%      RS.Yutm  
%      RS.Z % spatial filter avg elev NAVD88
%      RS.Sigma   % spatial filter standard deviation "roughness"

winsize=11; % 2d spatial window size in meters
minNR=5; % minimum number of data points in window to save a sigma value
RS.Mopnum=MopNumber;
RS.SpatialWindowSize=winsize;
RS.SpatialMinPoints=minNR;

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
[avgZ1,stdZ1]=runfilt2d(winsize,Z,minNR); %use winsize and monNP settings

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
[avgZ2,stdZ2]=runfilt2d(winsize,Z,minNR); % use winsize and monNP settings

%% blend grid data together.
%
R=nan(size(X));% initialize composite roughness grid as no data NaNs
D=nan(size(X));% initialize composite depth grid as no data NaNs
% 2009 has data, 2014 does not
idx=find(isnan(stdZ2(:)) & ~isnan(stdZ1(:))); 
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
load('MopTableUTM.mat','Mop');
N=XY2MopNumsV2(X(idx),Y(idx),Mop);
% retain only the x,y,z,roughness data in the input mop area
RS.Xutm=X(idx(N == MopNumber));
RS.Yutm=Y(idx(N == MopNumber));
RS.Z=D(idx(N == MopNumber));
RS.Sigma=R(idx(N == MopNumber));

end

