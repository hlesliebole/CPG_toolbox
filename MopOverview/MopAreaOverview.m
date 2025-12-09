% Mop area overview
clearvars
MopNumber=510;
fprintf('\nMop %i Mobile LiDAR Data Overview\n\n',MopNumber);

load(['M' num2str(MopNumber,'%5.5i') 'SA.mat']);
SA=SAcombineMops(505,518);
% survey indices of mobile LiDAR
%ldx=find((strcmp({SA.Source},'UTAir') | strcmp({SA.Source},'KMair')));
ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')));
%ldx=find((strcmp({SA.Source},'USACE')));
fprintf(' 1. %i total Mobile LiDAR surveys\n',numel(ldx));
fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));

%%

% Figure out total area within Mop area defined by standard Mop back beach
% and offshore points.
load MopTableUTM.mat

% 7 point (pt 1 = pt 7) UTM mop area boundary line starting a back beach 
% point and using mid-points between adjacent mop back/off points 
MopXutmBound(1)=Mop.BackXutm(MopNumber);
MopYutmBound(1)=Mop.BackYutm(MopNumber);
MopXutmBound(2)=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
MopYutmBound(2)=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
MopXutmBound(3)=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
MopYutmBound(3)=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
MopXutmBound(4)=mean(Mop.OffXutm(MopNumber));
MopYutmBound(4)=mean(Mop.OffYutm(MopNumber));
MopXutmBound(5)=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber-1)]);
MopYutmBound(5)=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber-1)]);
MopXutmBound(6)=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber-1)]);
MopYutmBound(6)=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber-1)]);
MopXutmBound(7)=MopXutmBound(1);
MopYutmBound(7)=MopYutmBound(1);

% [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
[MopLatBoundary,MopLonBoundary]=utm2deg(MopXutmBound,MopYutmBound,...
    repmat('11 S',[length(MopXutmBound) 1]));
% 
% x1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
% y1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
% x2=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
% y2=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
% x3=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
% y3=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
% x4=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
% y4=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);

%fill([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[.95 .95 .95],'edgecolor','none')
figure('position',[134   112   885   685]);
fill(MopLonBoundary,MopLatBoundary,[.95 .95 .95],'FaceAlpha',.2);
hold on;
plot([Mop.BackLon(MopNumber) Mop.OffLon(MopNumber)],...
    [Mop.BackLat(MopNumber) Mop.OffLat(MopNumber)],'m*-','linewidth',2);

plot([Mop.BackLon(MopNumber+1) Mop.OffLon(MopNumber+1)],...
    [Mop.BackLat(MopNumber+1) Mop.OffLat(MopNumber+1)],'y*-','linewidth',2);

plot([Mop.BackLon(MopNumber-1) Mop.OffLon(MopNumber-1)],...
    [Mop.BackLat(MopNumber-1) Mop.OffLat(MopNumber-1)],'y*-','linewidth',2);

plot_google_map('MapType', 'satellite','Alpha', 1);
set(gca,'FontSize',12)

set(gca,'xlim',[Mop.BackLon(MopNumber)-0.0015 Mop.BackLon(MopNumber)+0.0015]);
set(gca,'ylim',[Mop.BackLat(MopNumber)-0.0015 Mop.BackLat(MopNumber)+0.0015])
plot_google_map('MapType', 'satellite','Alpha', 1);

MobileLidarXutm=vertcat(SA(ldx).X);
MobileLidarYutm=vertcat(SA(ldx).Y);

% reduce to unique x,y points
[uxy,ia,ic]=unique([MobileLidarXutm,MobileLidarYutm],'rows');
MobileLidarXutm=MobileLidarXutm(ia);MobileLidarYutm=MobileLidarYutm(ia);

% create x,y boundary vector around all the lidar data points
% start at max X at min Y point
ymin=min(MobileLidarYutm);ymax=max(MobileLidarYutm);
np=0; % boudary point counter
% trace northward east boundary
for y=ymin:ymax
    iy=find(MobileLidarYutm == y); % LiDAR points with this y
    if ~isempty(iy)
     np=np+1;
     [xmax,imax]=max(MobileLidarXutm(iy)); % max x point on y line
     MobileLidarXutmBound(np)=MobileLidarXutm(iy(imax));
     MobileLidarYutmBound(np)=MobileLidarYutm(iy(imax));
    end
end
% if there is more than one lidar point at the ymax boundary start at ymax on
%   southwrd steps, otherwise start a ymax-1
if numel(iy) > 1;ystart=ymax;else;ystart=ymax-1;end

% trace southward west boundary
for y=ystart:-1:ymin
    iy=find(MobileLidarYutm == y); % LiDAR points with this y
    if ~isempty(iy)
     np=np+1;
     [xmin,imin]=min(MobileLidarXutm(iy)); % max x point on y line
     MobileLidarXutmBound(np)=MobileLidarXutm(iy(imin));
     MobileLidarYutmBound(np)=MobileLidarYutm(iy(imin));
    end
end

% if the last point is not equal to the first boundary point, close the
%  boundary loop
if (MobileLidarXutmBound(end) ~= MobileLidarXutmBound(1) || ...
        MobileLidarYutmBound(end) ~= MobileLidarYutmBound(1) )
    MobileLidarXutmBound(end+1)=MobileLidarXutmBound(1);
    MobileLidarYutmBound(end+1)=MobileLidarYutmBound(1);
end

% [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
[MobileLidarLatBound,MobileLidarLonBound]=...
    utm2deg(MobileLidarXutmBound,MobileLidarYutmBound,...
    repmat('11 S',[length(MobileLidarXutmBound) 1]));

% [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
[MobileLidarLat,MobileLidarLon]=...
    utm2deg(MobileLidarXutm,MobileLidarYutm,...
    repmat('11 S',[length(MobileLidarXutm) 1]));

% overlay on map
plot(MobileLidarLon,MobileLidarLat,'.');hold on;
plot(MobileLidarLonBound,MobileLidarLatBound,'-');

fill(MopLonBoundary,MopLatBoundary,[.95 .95 .95],'FaceAlpha',.2);
hold on;
plot([Mop.BackLon(MopNumber) Mop.OffLon(MopNumber)],...
    [Mop.BackLat(MopNumber) Mop.OffLat(MopNumber)],'m*-','linewidth',2);

plot([Mop.BackLon(MopNumber+1) Mop.OffLon(MopNumber+1)],...
    [Mop.BackLat(MopNumber+1) Mop.OffLat(MopNumber+1)],'y*-','linewidth',2);

plot([Mop.BackLon(MopNumber-1) Mop.OffLon(MopNumber-1)],...
    [Mop.BackLat(MopNumber-1) Mop.OffLat(MopNumber-1)],'y*-','linewidth',2);

%% Use accumarray to show area statistics

MobileLidarXutm=vertcat(SA(ldx).X);
MobileLidarYutm=vertcat(SA(ldx).Y);
MobileLidarZnavd88=vertcat(SA(ldx).Z);

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(MobileLidarXutm);
[uy, ~, yidx] = unique(MobileLidarYutm);

%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.')./zcount;
zmedian = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @median); % 
zstd = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @std); % 
zmin = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @min); % 
zmax = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @max); % 
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final shore box data vectors
xutm=ux(i);yutm=uy(j);
zavg=zavg(ii);
zcount=zcount(ii);
zmedian=zmedian(ii);
zmin=zmin(ii);
zmax=zmax(ii);
zstd=zstd(ii);

% convert to lat lon
[lat,lon]=utm2deg(xutm,yutm,repmat('11 S',[length(xutm) 1]));

%% plot survey coverage count
z=zcount;
% plot colored points
minz=min(z);%quantile(z,.05);
maxz=max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String='Survey Coverge Count';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
title('Mobile LiDAR Survey Coverage Count')
pause
%% plot survey median
z=zmedian;
% plot colored points
minz=0;%min(z);%quantile(z,.05);
maxz=3;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String='Elevation (m, NAVD88) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
title('Mobile LiDAR Surveys Median Elevations')
pause
%% plot max surface
z=zmax;
% plot colored points
minz=0;%min(z);%quantile(z,.05);
maxz=3;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String=' Elevation (m, NAVD88) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
title('Mobile LiDAR Surveys Global Maximum Elevations')
pause

%% plot min surface
z=zmin;
% plot colored points
minz=0;%min(z);%quantile(z,.05);
maxz=3;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String=' Elevation (m, NAVD88) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
title('Mobile LiDAR Surveys Global Minimum Elevations')
pause
%% plot std surface
z=zstd;
% plot colored points
minz=0;%min(z);%quantile(z,.05);
maxz=1.8;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String=' Standard Deviation (m, NAVD88) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
title('Mobile LiDAR Surveys Elevation Standard Deviations')


%figure
%[ScatterPlot,ColorBarPlot]=ColorScatter2d(xutm,yutm,zcount);


% k=boundary(MobileLidarXutm,MobileLidarYutm,1);
% figure;plot(MobileLidarXutm(k),MobileLidarYutm(k),'.-');

% identify all the 1m avg points in the Mop area that have historical data

% calculate temporally weighted mean elevation grid and global min-max 
% elevation grids

