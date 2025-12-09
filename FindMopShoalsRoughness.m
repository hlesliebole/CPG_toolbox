function RG=FindMopShoalsRoughness(MopNumber)

%%  returns "roughness" grid struct array RG of SHOALS spatial standard deviations 
%      RG.Mopnum
%      RG.SpatialWindowSize;
%      RG.SpatialMinPoints;
%      RG.Xutm 
%      RG.Yutm  
%      RG.Z % spatial filter avg elev
%      RG.Sigma   % spatial filter standard deviation "roughness"



winsize=11; % 2d spatial window size in meters
minNR=5; % minimum number of data points in window to save a sigma value
RG.Mopnum=MopNumber;
RG.SpatialWindowSize=winsize;
RG.SpatialMinPoints=minNR;

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

%% Process the 2009 and 2014 surveys separately and make gridded
%    results for both
%
% make 2d X,Y grid arrays
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);

%    2009
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

%    2014
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
[avgZ2,stdZ2]=runfilt2d(winsize,Z,minNR); %use winsize and monNP settings

%% blend together.
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


idx=find(~isnan(R(:))); % grid points with data
% assign them to nearest mops
load('MopTableUTM.mat','Mop');
N=XY2MopNumsV2(X(idx),Y(idx),Mop);
RG.Xutm=X(idx(N == MopNumber));
RG.Yutm=Y(idx(N == MopNumber));
RG.Yutm=Y(idx(N == MopNumber));
RG.Z=D(idx(N == MopNumber));
RG.Sigma=R(idx(N == MopNumber));

end

% 
% %% Code to calculate plot a composite of 2009 and 2014 SHOALS 1m spatial averaged survey 
% %   data for a reach of Mop points. Shoals surveys occurred over mutiple days in 2009 and 2014.
% 
% %% Combine 1m avg Mop SA struct arrays for a specified (starting,ending) Mop reach 
% 
% %CS=SAcombineMops(481,506);  % la jolla cove to the marine room
% %CS=SAcombineMops(665,715);  % cardiff SB to north of encinitas point/swamis
% %CS=SAcombineMops(658,686);  % cardiff SB 
% %CS=SAcombineMops(545,555);  % south torrey
% CS=SAcombineMops(505,512);  % lj shores south of pier
% %CS=SAcombineMops(623,633);  % del mar
% %CS=SAcombineMops(73,83);  % silver strand
% %CS=SAcombineMops(2,18);  % border field
% 
% 
% 
% %% find the USACE shoals surveys in the combined CS survey struct array
% 
% %% option for just the 2009 or 2014 survey dates
% iyr=2009;
% %iyr=2014;
% sidx=find(strcmp({CS.Source},'USACE') & year(datetime([CS.Datenum],'convertfrom','datenum')) == iyr);%%  
% 
% %% option to find all of the shoals surveys in both years 
% %sidx=find(strcmp({CS.Source},'USACE'));
% 
% %% get global x,y limits for the composite grid
% xmin=[];xmax=[];ymin=[];ymax=[];
% for n=sidx  % loop through shoals surveys
%  xmin=min([xmin [CS(n).X]']);xmax=max([xmax [CS(n).X]']);
%  ymin=min([ymin [CS(n).Y]']);ymax=max([ymax [CS(n).Y]']);
% end
% 
% %% make 2d X,Y grid arrays
% [X,Y]=meshgrid(xmin:xmax,ymin:ymax);
% % initialize composite elev grid as no data NaNs
% Z=nan(size(X)); 
% 
% %% progressively overlay Shoals 1m average data onto the grid 
% %  so any grid points with redundant info keep the most recent survey
% for n=sidx % loop through shoal survey dates
%     % 1d grid indices with data
%     gdx=sub2ind(size(X),round(CS(n).Y)-ymin+1,round(CS(n).X)-xmin+1); 
%     Z(gdx)=CS(n).Z; % overlay on grid
% end
% 
% %---------------------------------------------------
% % elev change signal 2 noise qc flags reef areas
% wsize=11; % odd number sq area for filter (must be >= 5)
% [avgZ,stdZ]=runfilt2d(wsize,Z,5);
% %stdZ(stdZ < 0.3)=NaN;
% ndx=find(~isnan(stdZ(:)) & Z(:) < 0.774);
% x=X(ndx);y=Y(ndx);z=-stdZ(ndx);az=avgZ(ndx);
% % reduce to points with positive MOP shoreline xshore locations
% %   to remove any estuary points
% [Nmop,xs]=UTM2MopxshoreX(x,y);
% x(xs < 0)=[];y(xs < 0)=[];z(xs < 0)=[];az(xs < 0)=[];
% 
% figure;ScatterPlotSubaqueousUTM(x,y,z,min(z),'2d');
% figure;plot(az,-z,'k.')
% 
% end
