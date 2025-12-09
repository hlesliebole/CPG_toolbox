% Uses linear scatter interpolation/extrapolation of valid points to fill 
%  in subarial survey area

%% ------------------------------------------------------------------
%% 1. remove all the data at spatial grid points that are in less than
%%   10% of the surveys (mostly low points on beach and water points)
MinFractionSurveysWithData=0.1;

MobileLidarXutm=vertcat(SA(ldx).X);
MobileLidarYutm=vertcat(SA(ldx).Y);
%MobileLidarZnavd88=vertcat(SA(ldx).Z);

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
%MobileLidarZnavd88(nbad)=[];



%% ------------------------------------------------------------------
%% calculate global subaerial boundary line around all the valid
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

% % [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
% [MobileLidarLatBound,MobileLidarLonBound]=...
%     utm2deg(MobileLidarXutmBound,MobileLidarYutmBound,...
%     repmat('11 S',[length(MobileLidarXutmBound) 1]));

% % [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
% [MobileLidarLat,MobileLidarLon]=...
%     utm2deg(MobileLidarXutm,MobileLidarYutm,...
%     repmat('11 S',[length(MobileLidarXutm) 1]));
%% ------------------------------------------------------------------
%%  define all 1m spatial res Xutm,Yutm points on and within the subaerial
%%   boundary line.  This is the fixed area for all the subaerial
%%   volume estimates at this mop.

%% ------------------------------------------------------------------
%%  Loop through individual surveys and interpolate-extrapolate in
%%  2D to fill the subaerial mop area. 


sn=CurrentSurveyNumber;
x=vertcat(SA(ldx(sn)).X);
y=vertcat(SA(ldx(sn)).Y);
z=vertcat(SA(ldx(sn)).Z);

% remove any QC flagged points
xbad=vertcat(QC(ldx(sn)).X);
ybad=vertcat(QC(ldx(sn)).Y);
zbad=vertcat(QC(ldx(sn)).Z);
fprintf(' Removing %i bad points flagged by QC\n',numel(xbad))
if numel(xbad) > 0
 nbad=find(ismember([x';y']',[xbad;ybad]','rows'));
 fprintf(' Matched %i bad points flagged by QC\n',numel(nbad))
 x(nbad)=[];y(nbad)=[];z(nbad)=[];
end

% MobileLidarXutm=x;
% MobileLidarYutm=y;
% MobileLidarZnavd88=z;

    % create scattered interpolant function with survey data
    F1 = scatteredInterpolant(x,y,z,'linear','linear');
    % interpolate-extrapolate to get elevations at all points within
    %  fixed subaerial mop area boundary
    Zarea = F1(MobileLidarXunq,MobileLidarYunq);

%% ------------------------------------------------------------------
%% calculate the mop subaerial area volume above specified contour elevation
%%

%     ivol=find(Zarea >= VolumeElevCutoff);%1.344);
%     ilow=find(Zarea < VolumeElevCutoff);
%     if ~isempty(ivol)
%         ns=ns+1;
%         Sdate(ns)=datetime(SA(s).Datenum,'ConvertFrom','datenum');
%         VolMsl(ns)=sum(Zarea(ivol)-VolumeElevCutoff);
%         VolN(ns)=numel(find(MobileLidarZnavd88(idx) >= VolumeElevCutoff))/numel(ivol);%numel(idx);
%         Bmsl(ns)=numel(ilow);
%     end
% 
% z=Zarea;

%% Make scatter plot of gridded data

if exist('nav','var');delete(nav);end
if exist('ax1','var');delete(ax1);end

x=MobileLidarXunq;
y=MobileLidarYunq;
z=Zarea;

ax1=axes('position',[.1 .1 .8 .8]);
% plot colored points
minz=min(z);%quantile(z,.05);
maxz=max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
if dview == '3'
scp=scatter3(x,y,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 10 1]);
zoom off
elseif dview == 'T'
scp=scatter(x,y,12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 10 1]);
zoom on
elseif dview == 'N'
scp=scatter(x,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 1 1]);
zoom on
elseif dview == 'E'
scp=scatter(y,z, 12, cm(ceil(zscaled(idx)),:), 'filled');
set(ax1,'dataaspectratio',[10 1 1],'xdir','reverse');
zoom on
end

cb=colorbar;cb.Label.String='Elevation';
set(ax1,'clim',[minz maxz],'fontsize',16);colormap(cm);
title(['MOP ' num2str(CurrentMopNumber) ' ' datestr(SA(ldx(sn)).Datenum) ' ' SA(ldx(sn)).Source] )
set(ax1,'color',[.5 .5 .5])
