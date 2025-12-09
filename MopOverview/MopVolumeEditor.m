% Mop area overview
clearvars
close all
MopNumber=635%54%655;%768;%654;
MinFractionSurveysWithData=0.05;
MinFractionalCoverageAboveElevCutoff=0.5;
VolumeElevCutoff=1.344;%0.774;%-0.058;%1.344;%0.774;%

fprintf('\nMop %i Mobile LiDAR Data Overview\n\n',MopNumber);

load(['M' num2str(MopNumber,'%5.5i') 'SA.mat']);

ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')));
fprintf(' 1. %i total Mobile LiDAR surveys\n',numel(ldx));
fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));
MobileLidarXutm=vertcat(SA(ldx).X);
MobileLidarYutm=vertcat(SA(ldx).Y);
MobileLidarZnavd88=vertcat(SA(ldx).Z);
% make a survey date to go with every survey point
for n=1:size(SA,2)
    SA(n).Dates=repmat(SA(n).Datenum,size(SA(n).X));
end
MobileLidarDatenum=vertcat(SA(ldx).Dates);
for n=1:size(SA,2)
    SA(n).Survnum=repmat(n,size(SA(n).X));
end
MobileLidarSurvnum=vertcat(SA(ldx).Survnum);

[uxy,ia,ic]=unique([MobileLidarXutm,MobileLidarYutm],'rows');
fprintf(' 4. %i unique LiDAR Xutm,Yutm points in area.\n',numel(ia));

%% ------------------------------------------------------------------
%% 1. remove all the data at spatial grid points that are in less than
%%   10% of the surveys (mostly low points on beach and water points)

% bin all survey data points by placing in unique x,y data array
[ux, ~, xidx] = unique(MobileLidarXutm);
[uy, ~, yidx] = unique(MobileLidarYutm);

% array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  

% reduce arrays to 1d vectors of x,y points with accumarray z data 
ii=zcount > 0; % 1d indices of valid data
[i,j]=find(zcount > 0); % 2d indices of valid data
% final data vectors
xutm=ux(i);yutm=uy(j); % utm points with data
zcount=zcount(ii);

% remove lidar grid points that have data in less that 10% of the surveys
ibad=find(zcount < MinFractionSurveysWithData*numel(ldx));

nbad=find(ismember([MobileLidarXutm';MobileLidarYutm']',[xutm(ibad)';yutm(ibad)']','rows'));
MobileLidarXutm(nbad)=[];
MobileLidarYutm(nbad)=[];
MobileLidarZnavd88(nbad)=[];
MobileLidarDatenum(nbad)=[];
MobileLidarSurvnum(nbad)=[];

%% ------------------------------------------------------------------
%% calculate global subaerial boundary line around all the remaining
%%  mop area lidar points

% reduce to unique x,y points
[uxy,ia,ic]=unique([MobileLidarXutm,MobileLidarYutm],'rows');
MobileLidarXunq=MobileLidarXutm(ia);MobileLidarYunq=MobileLidarYutm(ia);

% create x,y boundary vector around all the lidar data points
% start at max X at min Y point
ymin=min(MobileLidarYunq);ymax=max(MobileLidarYunq);
np=0; % boudary point counter
% trace northward east boundary
for y=ymin:ymax
    iy=find(MobileLidarYunq == y); % LiDAR points with this y
    if ~isempty(iy)
     np=np+1;
     [xmax,imax]=max(MobileLidarXunq(iy)); % max x point on y line
     MobileLidarXutmBound(np)=MobileLidarXunq(iy(imax));
     MobileLidarYutmBound(np)=MobileLidarYunq(iy(imax));
    end
end
% if there is more than one lidar point at the ymax boundary start at ymax on
%   southwrd steps, otherwise start a ymax-1
if numel(iy) > 1;ystart=ymax;else;ystart=ymax-1;end

% trace southward west boundary
for y=ystart:-1:ymin
    iy=find(MobileLidarYunq == y); % LiDAR points with this y
    if ~isempty(iy)
     np=np+1;
     [xmin,imin]=min(MobileLidarXunq(iy)); % max x point on y line
     MobileLidarXutmBound(np)=MobileLidarXunq(iy(imin));
     MobileLidarYutmBound(np)=MobileLidarYunq(iy(imin));
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
%% ------------------------------------------------------------------
%%  define all 1m spatial res Xutm,Yutm points on and within the subaerial
%%   boundary line.  This is the fixed area for all the subaerial
%%   volume estimates at this mop.


%%

%%

figure('position',[5          72        1426         717]);

%%

Mhead=uicontrol(MainFig,'style','text','position',[10 660 70 30],...
    'string','MOP #','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',22);

Medit=uicontrol(MainFig,'style','edit','position',[85 660 100 30],...
    'string',CurrentMopNumber,'Value',1,'fontsize',22,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'CurrentMopNumber=str2num(Medit.String);delete(Mhead);delete(Medit);GetMopMobileLidarSA');



ax1=axes('position',[.05 .55 .25 .4]);
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

BeachAreaLon=mean(MobileLidarLonBound);BeachAreaLat=mean(MobileLidarLatBound);
set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
set(gca,'xtick',[],'ytick',[])
plot_google_map('MapType', 'satellite','Alpha', 1);
title({['MOP ' num2str(MopNumber) ' : Mobile LiDAR Subaerial Boundary'],....
    ['Minimum Frequency of Coverage: ' num2str(100*(MinFractionSurveysWithData),'%4.1f') '% of Surveys']},...
    'fontsize',14)

%% --------------------------------------------------------------------
%% Use accumarray to get area statistics

% MobileLidarXutm=vertcat(SA(ldx).X);
% MobileLidarYutm=vertcat(SA(ldx).Y);
% MobileLidarZnavd88=vertcat(SA(ldx).Z);
% reduce to MSL zone
% imsl=find(MobileLidarZnavd88 > 0.75 & MobileLidarZnavd88 < 1.0);
% MobileLidarXutm=MobileLidarXutm(imsl);
% MobileLidarYutm=MobileLidarYutm(imsl);
% MobileLidarZnavd88=MobileLidarZnavd88(imsl);

% bin all survey data points by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(MobileLidarXutm);
[uy, ~, yidx] = unique(MobileLidarYutm);

% array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.')./zcount;
zmedian = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @median); % 
zmode = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @mode); % 
zstd = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @std); % 
zmin = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @min); % 
zmax = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @max); % 
zq75 = accumarray([xidx(:), yidx(:)], MobileLidarZnavd88.',[], @(x) quantile(x,.75));
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)],MobileLidarZnavd88.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with accumarray z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final data vectors
xutm=ux(i);yutm=uy(j); % utm points with data
zavg=zavg(ii);
zcount=zcount(ii);
zmedian=zmedian(ii);
zmin=zmin(ii);
zmax=zmax(ii);
zstd=zstd(ii);
zmode=zmode(ii);
zq75=zq75(ii);
%%

ax2=axes('position',[.05 .05 .25 .4]);
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


% convert to lat lon
[lat,lon]=utm2deg(xutm,yutm,repmat('11 S',[length(xutm) 1]));

% 
%% plot min surface
z=zmin;
% plot colored points
minz=-0.5;%min(z);%quantile(z,.05);
maxz=5;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String=' Elevation (m, NAVD88) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
set(gca,'xtick',[],'ytick',[])
plot_google_map('MapType', 'satellite','Alpha', 1);

title('Mobile LiDAR Global Minimum Elevation Surface','fontsize',14)

%% ------------------------------------------------------------------
%%  Loop through individual surveys and interpolate-extrapolate in
%%  2D to fill the subaerial mop area. 


ns=0;
for s=ldx
    % find indices of data points for this survey
    %idx=find(SA(s).Datenum == MobileLidarDatenum);

    % if more than one lidar survey on same day, use the last one
    idx=find(s == MobileLidarSurvnum);
    %fprintf('%i\n',s)
    %numel(idx)
    % create scattered interpolant function with survey data
    if numel(idx) > 3
    F1 = scatteredInterpolant(MobileLidarXutm(idx),MobileLidarYutm(idx),...
        MobileLidarZnavd88(idx),'linear','linear');
    % interpolate-extrapolate to get elevations at all points within
    %  fixed subaerial mop area boundary
    Zarea = F1(MobileLidarXunq,MobileLidarYunq);
    else
        Zarea=[];
    end
%% ------------------------------------------------------------------
%% calculate the mop subaerial area volume above specified contour elevation
%%

if numel(Zarea) > 0
    ivol=find(Zarea >= VolumeElevCutoff);%1.344);
    ilow=find(Zarea < VolumeElevCutoff);
    itot=find(Zarea > zmin);
    if ~isempty(ivol)
        ns=ns+1;
        Sdate(ns)=datetime(SA(s).Datenum,'ConvertFrom','datenum');
        VolCutElev(ns)=sum(Zarea(ivol)-VolumeElevCutoff);
        VolMin(ns)=sum(Zarea(ivol)-zmin(ivol));
        VolN(ns)=numel(find(MobileLidarZnavd88(idx) >= VolumeElevCutoff))/numel(ivol);%numel(idx);
        Bmsl(ns)=numel(ilow);
        VolTot(ns)=sum(Zarea(itot)-zmin(itot));
    end
end

end

% remove low lidar coverage cases
VolCutElev(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;
VolMin(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;
VolTot(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;


ax3=axes('position',[.4 .55 .55 .4]);
p1=plot(Sdate,VolCutElev/1000,'k*-','linewidth',2);hold on
p2=plot(Sdate,VolMin/1000,'r*-','linewidth',2);hold on
p3=plot(Sdate,VolTot/1000,'g*-','linewidth',2);hold on
set(gca,'fontsize',16);grid on;
ylabel('MOP Area Volume (K m^{3})')
legend([p1 p2 p3],['Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88'],...
    ['Vol above Minimum Elev Surface, Landward of ' num2str(VolumeElevCutoff,'%5.3fm') ' contour'],...
    'Total Vol above Minimum Surface','location','north')
yl=get(gca,'ylim');set(gca,'ylim',[yl(1) yl(1)+(yl(2)-yl(1))*1.2]);
title(['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'])

%%
ax4=axes('position',[.4 .06 .55 .38]);
sy=year(Sdate);sm=month(Sdate);
hold on;
for y=min(sy):max(sy)
    if y == max(sy)
     plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'m*-',...
     'linewidth',3,'displayname',num2str(y));hold on
    else       
     plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'*-',...
     'linewidth',2,'displayname',num2str(y));hold on
    end
end
set(gca,'xlim',[1 13],'xtick',1:12,'fontsize',16);grid on;legend('location','eastoutside')
xlabel('Month of Year');ylabel('MOP Area Volume (K m^{3})')
title({['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'],...
    [' Seasonal view, Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88']})

