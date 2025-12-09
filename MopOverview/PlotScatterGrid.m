
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

%figure('position',[168   208   782   548]);
ns=0;
%for s=ldx
    s=ldx(CurrentSurveyNumber);
    % find indices of data points for this survey
    %idx=find(SA(s).Datenum == MobileLidarDatenum);
    idx=find(s == MobileLidarDatenum);
    %fprintf('%i\n',s)
    %numel(idx)
    % create scattered interpolant function with survey data
    F1 = scatteredInterpolant(MobileLidarXutm(idx),MobileLidarYutm(idx),...
        MobileLidarZnavd88(idx),'linear','linear');
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

ScatterPlotUTM(MobileLidarXunq,MobileLidarYunq,Zarea,'3d');


%end