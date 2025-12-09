%% Builds  the struct array TransectVols for a range of Mops
clearvars

MinFractionSurveysWithData=0.05;
MinFractionalCoverageAboveElevCutoff=0.5;
VolumeElevCutoff=0.774;%-0.058;%1.344;%0.774;%

Nmop=0;
for MopNumber=2:167%519:769%20:70%519:769
    CurrentMopNumber=MopNumber;
Nmop=Nmop+1;

% load the SA.mat file
SAmatfile=['M' num2str(CurrentMopNumber,'%5.5i') 'SA.mat'];
fprintf('Loading SA struct array from %s\n',SAmatfile)
load(SAmatfile);

% load the QC mat file if it exists
QCmatfile=['M' num2str(CurrentMopNumber,'%5.5i') 'QC.mat'];
% get QC struct array if exist, otherwise make one to mirror SA
%  but without any matching x,y,z qc removal points
fprintf('Loading QC struct array from %s\n',QCmatfile)
if exist(QCmatfile,'file')
    load(QCmatfile);
    % syncronize SA and QC survey dates, to catch source data file changes
    %   synced QC struct array is SQC
    SQC=SA; % start with mirror of SA
    for n=1:size(SA,2) 
        % match SA & QC reefbreak original lidar data file name and file modification date
        idx=find(strcmp({QC.File},SA(n).File) & [QC.FileDatenum] == SA(n).FileDatenum);
        if numel(idx) > 0
        SQC(n).X=QC(idx).X;SQC(n).Y=QC(idx).Y; SQC(n).Z=QC(idx).Z;
        else % if no direct match, no QC changes
        SQC(n).X=[]; SQC(n).Y=[]; SQC(n).Z=[];
        end
    end
    
else
    fprintf('QC file: %s\n Not found. Making empty synced QC struct array.\n',QCmatfile)
    SQC=SA;
    for n=1:size(SA,2)
        SQC(n).X=[]; SQC(n).Y=[]; SQC(n).Z=[];
    end
end

fprintf('\nMop %i Mobile LiDAR Data Overview\n\n',CurrentMopNumber);

%ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')));
ldx=find( ( strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR') | ...
    strcmp({SA.Source},'TrkMR') | strcmp({SA.Source},'UTAir') |...
    strcmp({SA.Source},'KMair')  | strcmp({SA.Source},'USGS')) );
% ldx=find( ( strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR') | ...
%     strcmp({SA.Source},'TrkMR') ) );
fprintf(' 1. %i total Mobile LiDAR surveys\n',numel(ldx));
fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));


%% Remove any QC flagged data points from the SA struct array
nqc=0; % qc point counter
for n=1:size(SQC,2) % loop through all surveys
    if numel(SQC(n).X) > 0 % check if a survey has bad data
        x=vertcat(SA(n).X);
        y=vertcat(SA(n).Y);
         %z=vertcat(SQC(n).Z);
        xbad=vertcat(SQC(n).X);
        ybad=vertcat(SQC(n).Y);
         %zbad=vertcat(SQC(n).Z);
         %fprintf(' Removing %i bad points flagged by QC\n',numel(xbad))
        % find the SA indices of the bad QC points and remove from SA
        nbad=find(ismember([x';y']',[xbad;ybad]','rows'));
        nqc=nqc+numel(nbad); % add to qc point total
          %fprintf(' Matched %i bad points flagged by QC\n',numel(nbad))
        SA(n).X(nbad)=[];
        SA(n).Y(nbad)=[];
        SA(n).Z(nbad)=[];
    end
end

fprintf(' %i SA data points removed using the synced QC struct array.\n',nqc)

MobileLidarXutm=vertcat(SA(ldx).X);
MobileLidarYutm=vertcat(SA(ldx).Y);
MobileLidarZnavd88=vertcat(SA(ldx).Z);

% make a survey date and survey number to go with every survey point
for n=1:size(SA,2)
    SA(n).Dates=repmat(SA(n).Datenum,size(SA(n).X));
end
MobileLidarDatenum=vertcat(SA(ldx).Dates);
for n=1:size(SA,2)
    SA(n).Survnum=repmat(n,size(SA(n).X));
end
MobileLidarSurvnum=vertcat(SA(ldx).Survnum);

%% special case for Mop 638 at dog beach. only keep points west of
%  easting 474705
if CurrentMopNumber == 638
    MobileLidarZnavd88=MobileLidarZnavd88(MobileLidarXutm < 474705);
    MobileLidarYutm=MobileLidarYutm(MobileLidarXutm < 474705);
    MobileLidarDatenum=MobileLidarDatenum(MobileLidarXutm < 474705);
    MobileLidarSurvnum=MobileLidarSurvnum(MobileLidarXutm < 474705);
    MobileLidarXutm=MobileLidarXutm(MobileLidarXutm < 474705);
end

%% special case for Mop 635 north del mar. only keep points west of
%  easting 474800
if CurrentMopNumber == 635
    iout=find(MobileLidarXutm > 474800);
    
    MobileLidarZnavd88(iout)=[];
    MobileLidarYutm(iout)=[];
    MobileLidarDatenum(iout)=[];
    MobileLidarSurvnum(iout)=[];
    MobileLidarXutm(iout)=[];
end


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
MobileLidarXutmBound=[];
MobileLidarYutmBound=[];

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
%vax1=axes('position',[.05 .55 .25 .35]);
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
%% --------------
%% ---------------
%% - create vectors of all the unique Xutm and Yutm 1m grid points within the 
%     Mop subarial boundary. Some of points might not have data may not
%     have past data, so create a grid of x,y values based on the min-max
%     range of the boundary points and use the inpolygon function to 
%     reduce to the unique points within the boundary

[Xg,Yg]=meshgrid(min(MobileLidarXutmBound):max(MobileLidarXutmBound),...
    min(MobileLidarYutmBound):max(MobileLidarYutmBound));
%figure;plot(MobileLidarXutmBound,MobileLidarYutmBound,'k-');
%hold on;
%plot(Xg(:),Yg(:),'r.')
idx=find(inpolygon(Xg(:),Yg(:),MobileLidarXutmBound,MobileLidarYutmBound));
%plot(Xg(idx),Yg(idx),'c.')
% redefine unque x,y to based on inpolygon results
MobileLidarXunq=Xg(idx);
MobileLidarYunq=Yg(idx);

%% - now relate the unique x,y points in the subaerial zone to their nearest
%     1m xshore resolution point on the MOP transect

 
%   Get the main transect as UTM xt,yt points and subtransect lines as 
%   UTM xst,yst points ~1m apart.  This also sets the fixed 1m spaced X1D transect 
%   points starting -100m behind the mop back beach point out to the offshore point
NumSubTrans=1;
[X1D,xt,yt,xst,yst]=GetTransectLines(Mop,CurrentMopNumber,NumSubTrans,[-100 0]);

clear Zdist Znear
    
    % find the distances of all data apoints to the main transect line points  
     [dp,NearIdx]=pdist2([yt',xt'],...
         [double(MobileLidarYunq),double(MobileLidarXunq)],...
        'euclidean','smallest',1);
    
    % define X based on the nearest transect line point
    %   row=nearest subtransect number; col = xshore distance indice on
    %   the nearest subtransect
     [row,col] = ind2sub(size(xt),NearIdx); 
    
    % xshore distance (1m xshore resolution) along the Mop transect for each survey point
     X=X1D(col);
    
    % For the xshore distances (1m xshore resolution) with data (Xuniq), 
    %  find the nearest survey point elevation for each (Znear) and
    %  how far away from the line it was (Zdist).
    
     Xuniq=unique(X);
    
     n=0;
     for x=Xuniq
        n=n+1; % xshore point counter
        Transect.X1D=x;
        idx=find(X == x);
        % Transect.Xutm=MobileLidarXunq(idx);
        % Transect.Yutm=MobileLidarYunq(idx);
        % Transect.UniqueXYindexes=idx;
        % hold on;plot(MobileLidarXunq(idx),MobileLidarYunq(idx),'.');
        [Zdist(n),imin]=min(dp(idx));
        Znear(n)=1; % assign an elevation of 1 as valid transect points  
     end
    
     % set nearest point elevations with Zdist > 1m to NaNs
     Znear(Zdist > 1.5)=NaN;
    
    % seed the regular spaced profile vectors (z1d) with the 
    %  valid nearest point elevations (both valid and NaNs)
     ndx=find(ismember(X1D,Xuniq));
     z1d=X1D*NaN;z1d(ndx)=Znear;

     % find the X1D back beach boundary for this subtransect
     %  as the first valid transect point with lidar data within 1m
     idxBack=find(~isnan(z1d), 1 );
     Xback=X1D(idxBack);

     % assign unique grid points for transect points behind the 
     % lidar back beach point to the back point.
     n=0;
     A1D=X1D*0;
     % loop through transect points with subaerial mop area grid points
     %  associated with them.
     for x=Xuniq     
        if x >= Xback
         n=n+1; % xshore point counter
         MopTransectArea(n).X1D=x;
          idxX=find(X1D == x);
          MopTransectArea(n).Xutm=xt(idxX);
          MopTransectArea(n).Yutm=yt(idxX);
          if x == Xback
            idx=find(X <= x);
          else
            idx=find(X == x);
          end
          % number of Mop Area 1m grid points associated with a 1m transect
          %  line point
          A1D(idxX)=numel(idx); 
         MopTransectArea(n).Xutms=MobileLidarXunq(idx);
         MopTransectArea(n).Yutms=MobileLidarYunq(idx);
         MopTransectArea(n).UniqueXYindexes=idx;
         hold on;plot(MobileLidarXunq(idx),MobileLidarYunq(idx),'.');
         plot(MopTransectArea(n).Xutm,MopTransectArea(n).Yutm,'k.')
        end
     end
     
%% --------------------------------------------------------------------
%  Calculate Volumes for every transect calculated from SA

%load(SAmatfile);
NumSubTrans=1;
XgapTol=5;
YdistTol=50;
[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

% remove any elevations seaward of the min profile elevations
    %  (lidar swash filter)
     Zf=Z1Dmedian;
     [zmin,imin]=min(Zf');
    for n=1:size(Zf,1)
      if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
      end
    end


% remove any additional outliers
TF=isoutlier(Zf,"mean");
Zf(TF == 1)=NaN;

ZminGlobal=min(Zf);

ia1d=find(A1D > 0);
Sdate=datetime([],[],[]);
VolCutElev=[];
VolMin=[];
VolTot=[];
ns=0;
for s=1:size(SA,2)

    igood=find(~isnan(Zf(s,:)));
    if numel(igood) > 10
    z1d=interp1(X1Dmop(igood),Zf(s,igood),X1Dmop(ia1d),'linear','extrap');
    z1dmin=interp1(X1Dmop(igood),ZminGlobal(igood),X1Dmop(ia1d),'linear','extrap');
    
    ivol=find(z1d >= VolumeElevCutoff);%1.344);
    itot=find(z1d > z1dmin);
    if ~isempty(ivol)
        ns=ns+1;
        Sdate(ns)=datetime(SA(s).Datenum,'ConvertFrom','datenum');
        VolCutElev(ns)=sum((z1d(ivol)-VolumeElevCutoff).*A1D(ia1d(ivol)));
        VolMin(ns)=sum((z1d(ivol)-z1dmin(ivol)).*A1D(ia1d(ivol)));
        VolTot(ns)=sum((z1d(itot)-z1dmin(itot)).*A1D(ia1d(itot)));
    end
    end
end

%% -----------------
%% ----------------

%vax3=axes('position',[.4 .55 .5 .35]);
% plot using datenum so ginput function works on it
Sdatenum=datenum(Sdate);
%numel(find(~isnan(VolTot)))
fprintf('\n%i valid volume estimates.\n\n',numel(find(~isnan(VolTot))));

TransectVols(Nmop).Mop=CurrentMopNumber;
TransectVols(Nmop).CutElev=VolumeElevCutoff;
TransectVols(Nmop).TotAreaMinElev=[];
TransectVols(Nmop).Vdatetimes=datetime([],[],[]);
TransectVols(Nmop).Vdatenums=[];
TransectVols(Nmop).VolCutElev=[];
TransectVols(Nmop).VolMin=[];
TransectVols(Nmop).VolTot=[];


if numel(Sdatenum) > 0

TransectVols(Nmop).Vdatetimes=Sdate;
TransectVols(Nmop).TotAreaMinElev=min(zmin(:));
TransectVols(Nmop).Vdatetimes=Sdate;
TransectVols(Nmop).Vdatenums=Sdatenum;
TransectVols(Nmop).VolCutElev=VolCutElev;
TransectVols(Nmop).VolMin=VolMin;
TransectVols(Nmop).VolTot=VolTot;

end

end

save MopTransectrVolumesMSLIB.mat TransectVols 

%%
%MobileTransectVols=vertcat(TransectVols.VolMin);


% p3=plot(Sdatenum,VolTot/1000,'g*-','linewidth',2);hold on
% p2=plot(Sdatenum,VolMin/1000,'r*-','linewidth',2);hold on
% p1=plot(Sdatenum,VolCutElev/1000,'k*-','linewidth',2);hold on
% 
% set(gca,'xtick',datenum(year(Sdate(1)):(year(Sdate(end))+1),1,1))
% datetick('x','yyyy','keepticks')
% set(gca,'fontsize',14);grid on;
% ylabel('MOP Area Volume (K m^{3})')
% legend([p1 p2 p3],['Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88'],...
%     ['Vol above Minimum Elev Surface, Landward of ' num2str(VolumeElevCutoff,'%5.3fm') ' contour'],...
%     'Total Vol above Minimum Surface','location','north')
% yl=get(gca,'ylim');set(gca,'ylim',[yl(1) yl(1)+(yl(2)-yl(1))*1.2]);
% title(['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'])
% 
% zoom on;
% option to edit a survey
% Vedit=uicontrol(gcf,'Style','Pushbutton','Position',[580 650 150 25],...
% 'Callback','GetVolSurveyNumber;EditSA','String','Edit a Survey'); 
% set(Vedit,'backgroundcolor','yellow')


%%
%vax4=axes('position',[.4 .08 .55 .35]);
% sy=year(Sdate);sm=month(Sdate);
% hold on;
% for y=min(sy):max(sy)
%     if y == max(sy)
%      plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'m*-',...
%      'linewidth',3,'displayname',num2str(y));hold on
%     else       
%      plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'*-',...
%      'linewidth',2,'displayname',num2str(y));hold on
%     end
% end
% set(gca,'xlim',[1 13],'xtick',1:12,'fontsize',14);grid on;legend('location','eastoutside')
% xlabel('Month of Year');ylabel('MOP Area Volume (K m^{3})')
% title({['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'],...
%     [' Seasonal view, Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88']})

