%  Example Code to blend various gridded survey data from
%  specific surveys/dates, in a specific order from older to newer, 
%  to make a best guess data grid.
% 
%  Final grid is in UTM coords with 1m spatial resolution
%   with navd88 elevations (MSL = 0.774 navd88)
%
%   Does not include SIO surveys of areas landward of the beach
%

%% add path to cpg mop files on reefbreak1 here
%  eg. on a mac with reefbreak1 mounted
addpath /volumes/group/mops
addpath /volumes/group/mops/toolbox

%% grid settings
% Define alongshore Mop range of grid area.
Mop1=560;
Mop2=620;
% define desired min depth (m) of offshore boundary 
DepthMin=-50;
% define approximate distance inland the grid should go (m)
DistInland=500;

%% load MOP transect info to get alongshore grid range
load('MopTableUTM.mat','Mop');

%% Work in lat lon to find the NOAA grid data of interest
%% set utm y range based on mop locations
ymin=ceil(min([Mop.BackYutm(Mop1) Mop.OffYutm(Mop1)]));
ymax=floor(max([Mop.BackYutm(Mop2) Mop.OffYutm(Mop2)]));

%% find the max east utm location of mop back beach points
%   to define the inland boundary
xmax=floor(max(Mop.BackXutm(Mop1:Mop2)))+DistInland;

%% -----------------------------------------------------------------
%    STEP 1: Make a base grid using NOAA PMEL DEM
% ------------------------------------------------------------------

%% load the NOAA PMEL MHW 1/3 arc sec data for the san diego region
%   (san_diego_13_mhw_2012.nc is in mops/toolbox)
fprintf('Loading san_diego_13_mhw_2012.nc basemap...\n')
xn=ncread('san_diego_13_mhw_2012.nc','lon');
yn=ncread('san_diego_13_mhw_2012.nc','lat');
zn=ncread('san_diego_13_mhw_2012.nc','Band1');
%  z is relative to MHW; subtract -1.344m to adjust to NAVD88
zn=zn-1.344;
% 
% % convert xn,yn lon,lat vectors into 2d indices

fprintf('Converting to UTM...\n')
% convert the 1d lat,lon grid indices to utm
%  using the mean latitude of the mop reach
AvgLat=mean(Mop.BackLat(Mop1:Mop2));
[xnUTM,y,utmzone]=deg2utm(AvgLat*ones(numel(xn),1),xn(:));
[x,ynUTM,utmzone]=deg2utm(yn,xn(1)*ones(numel(yn),1));

% There is a N and E offset of the PMEL data relative to the SIO
%  data. Shift ynUTM values south 275m and west 25m
ynUTM=ynUTM-275;
xnUTM=xnUTM-25;

%% trim NOAA grid info using the y limits and max x limit
%   keep and edge of extra data for regridding into 1m
%    spatial resolution within the min-max limits
edge=200; % 200m edge
idx=find(xnUTM <= xmax+edge);
idy=find(ynUTM >= ymin-edge & ynUTM <= ymax+edge);
xnUTM=xnUTM(idx);ynUTM=ynUTM(idy);zn=zn(idx,idy);

%% find the min east utm location in the NOAA grid data that 
%   is deeper than DepthMin setting
[ix,iy]=find(zn > DepthMin);
xmin=ceil(xnUTM(min(ix)));
idx=find(xnUTM >= xmin-edge);
xnUTM=xnUTM(idx); % trim x dimension
zn=zn(idx,:); % trim elevation grid

%% Now regrid the NOAA elevations to 1m spatial resolution
%    to form the base grid for blending
[yu,xu]=meshgrid(ynUTM,xnUTM); % noaa grid 2d indices
[y,x]=meshgrid(ymin:ymax,xmin:xmax); % desired 1m 2d indices
z=griddata(xu,yu,zn,x,y); % regrid with 1m spatial resolution

% make a figure of 1m spatial res base grid
figure('position',[ 17   241   537   543]);
imagesc(xmin:xmax,ymin:ymax,(z-.774)'); % msl plot
demcmap(z);cb=colorbar;cb.Label.String='Elev. (m, MSL)';
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('NOAA PMEL Base Grid');xlabel('E UTM');ylabel('N UTM');

%% -----------------------------------------------------------------
%    STEP 2: Make a combined gridded survey file struct array
%            for the defined Mop range
% ------------------------------------------------------------------
fprintf('Loading gridded Mop data...\n')
% get a extra mop on either side of grid boundaries 
SG=SGcombineMops(Mop1-1,Mop2+1);

% SG(N) contains 1m UTM spatial resolution gridded survey values for
%  each of the N historical surveys that had data between Mop1 
%  and Mop2.  Only points with valid gridded data are saved in
%  the SG(N).X .Y and .Z fields

%% -----------------------------------------------------------------
%    STEP 3: First replace noaa base grid values with USACE shoals
%            topobathy survey gridded data as this data goes
%            deeper than the jetski. 
% ------------------------------------------------------------------
fprintf('Overlaying gridded USACE SHOALS topobathy liDAR data...\n')
% find USACE Shoals survey indices
ShoalsIdx=find(strcmp({SG.Source},'USACE'));

% They are ordered oldest to newest so loop through
%  gridded shoals survey data and replace any 
%  base grid points with valid shoals gridded points. 

for sn=ShoalsIdx(2:end)
  % 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  % gridded shoals data cab be weird around the surfzone
  %  and is mostly wanted for deeper depths so just use 
  %  data deeper than -5m
  
  idx=find( xs > xmin & xs < xmax & ...
      ys > ymin & ys < ymax & ...
      zs < -5);
  
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 

end

% make a figure of 1m spatial res base grid
figure('position',[ 27   231   537   543]);
imagesc(xmin:xmax,ymin:ymax,(z-.774)');
demcmap(z);cb=colorbar;cb.Label.String='Elev. (m, MSL)';
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals');xlabel('E UTM');ylabel('N UTM');

%% -----------------------------------------------------------------
%    STEP 4: Now overlay jumbo surveys up to Jan 24 2023 
% ------------------------------------------------------------------
fprintf('Overlaying gridded SIO jumbo survey data...\n')
JumboIdx=find(contains({SG.File},'umbo') & ...
    [SG.Datenum] > datenum(2019,1,1) & ...
    [SG.Datenum] <= datenum(2023,1,24));

for sn=JumboIdx
    
  %[XGutm,YGutm,Z]=UTMalongshoreGrid(582,172);
% 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  
  % Mostly want subaqueous and swash zone data so keep everything
  %  below +1m navd88
  idx=find( xs > xmin & xs < xmax & ...
      ys > ymin & ys < ymax & ...
      zs < 1);
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 
end

figure('position',[ 37   221   537   543]);
imagesc(xmin:xmax,ymin:ymax,(z-.774)'); % msl plot
demcmap(z);cb=colorbar;cb.Label.String='Elev. (m, MSL)';
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals + Jumbos');xlabel('E UTM');ylabel('N UTM');

%% -----------------------------------------------------------------
%    STEP 5: Now overlay Jan 17 to  Feb 2 truck surveys (2 feb includes blacks) 
% ------------------------------------------------------------------
fprintf('Overlaying SIO truck LiDAR survey data...\n')
TrkIdx=find(strcmp({SG.Source},'Trk') & ...
    [SG.Datenum] > datenum(2023,1,16) & ...
    [SG.Datenum] <= datenum(2023,2,2));

for sn=TrkIdx
% 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  
  idx=find( xs > xmin & xs < xmax & ...
      ys > ymin & ys < ymax);
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 
end

figure('position',[ 37   221   537   543]);
imagesc(xmin:xmax,ymin:ymax,(z-.774)'); % msl plot
demcmap(z);cb=colorbar;cb.Label.String='Elev. (m, MSL)';
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals + Jumbos + Trk');xlabel('E UTM');ylabel('N UTM');

%% -----------------------------------------------------------------
%    STEP 6: Now overlay Pinc jetski surveys  
% ------------------------------------------------------------------
fprintf('Overlaying gridded PiNC jetski survey data...\n')
PincIdx=find(contains({SG.File},'PiNC')); 

for sn=PincIdx
% 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  
  idx=find( xs > xmin & xs < xmax & ...
      ys > ymin & ys < ymax);
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 
end

figure('position',[ 37   221   537   543]);
imagesc(xmin:xmax,ymin:ymax,(z-.774)'); % msl plot
demcmap(z);cb=colorbar;cb.Label.String='Elev. (m, MSL)';
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals + Jumbos + Trk + PiNC Jetskis');xlabel('E UTM');ylabel('N UTM');