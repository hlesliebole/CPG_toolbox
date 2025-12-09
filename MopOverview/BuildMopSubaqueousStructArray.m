%% Builds  the struct array LidarVols for a range of Mops
close all
clearvars

% MinFractionSurveysWithData=0.05;
% MinFractionalCoverageAboveElevCutoff=0.5;
VolumeElevCutoff=0.774;%-0.058;%1.344;%0.774;%
VolumeElevCutoff=1.344%2.5;

Nmop=0;
for MopNumber=584%519:684%519:769
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

fprintf('\nMop %i Jumbo Data Overview\n\n',CurrentMopNumber);

ldx=find(contains({SA.File},'umbo'));

fprintf(' 1. %i total Jumbo surveys\n',numel(ldx));
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

JumboXutm=vertcat(SA(ldx).X);
JumboYutm=vertcat(SA(ldx).Y);
JumboZnavd88=vertcat(SA(ldx).Z);

% make a survey date and survey number to go with every survey point
for n=1:size(SA,2)
    SA(n).Dates=repmat(SA(n).Datenum,size(SA(n).X));
end
JumboDatenum=vertcat(SA(ldx).Dates);
for n=1:size(SA,2)
    SA(n).Survnum=repmat(n,size(SA(n).X));
end
JumboSurvnum=vertcat(SA(ldx).Survnum);

%% ------------------------------------------------------------------
%%  define all 1m spatial res Xutm,Yutm points on and within the subaerial
%%   boundary line.  This is the fixed area for all the subaerial
%%   volume estimates at this mop.
%%
%figure('position',[100 100 1000 800]);
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

% fill(MopLonBoundary,MopLatBoundary,[.95 .95 .95],'FaceAlpha',.2);
% hold on;
% plot([Mop.BackLon(MopNumber) Mop.OffLon(MopNumber)],...
%     [Mop.BackLat(MopNumber) Mop.OffLat(MopNumber)],'m*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber+1) Mop.OffLon(MopNumber+1)],...
%     [Mop.BackLat(MopNumber+1) Mop.OffLat(MopNumber+1)],'y*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber-1) Mop.OffLon(MopNumber-1)],...
%     [Mop.BackLat(MopNumber-1) Mop.OffLat(MopNumber-1)],'y*-','linewidth',2);
% 
% %plot_google_map('MapType', 'satellite','Alpha', 1);
% set(gca,'FontSize',12)
%%
% % overlay on map
% plot(JumboLon,JumboLat,'.');hold on;
% plot(JumboLonBound,JumboLatBound,'-');
% 
% fill(MopLonBoundary,MopLatBoundary,[.95 .95 .95],'FaceAlpha',.2);
% hold on;
% plot([Mop.BackLon(MopNumber) Mop.OffLon(MopNumber)],...
%     [Mop.BackLat(MopNumber) Mop.OffLat(MopNumber)],'m*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber+1) Mop.OffLon(MopNumber+1)],...
%     [Mop.BackLat(MopNumber+1) Mop.OffLat(MopNumber+1)],'y*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber-1) Mop.OffLon(MopNumber-1)],...
%     [Mop.BackLat(MopNumber-1) Mop.OffLat(MopNumber-1)],'y*-','linewidth',2);
% 
% BeachAreaLon=mean(JumboLonBound);BeachAreaLat=mean(JumboLatBound);
% set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
% set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
% set(gca,'xtick',[],'ytick',[])
% plot_google_map('MapType', 'satellite','Alpha', 1);
% title({['MOP ' num2str(MopNumber) ' : Mobile LiDAR Subaerial Boundary'],....
%     ['Minimum Frequency of Coverage: ' num2str(100*(MinFractionSurveysWithData),'%4.1f') '% of Surveys']},...
%     'fontsize',14)

%% - create vectors of all the unique Xutm and Yutm 1m grid points within the 
%     Mop subarial boundary. Some of points might not have data may not
%     have past data, so create a grid of x,y values based on the min-max
%     range of the boundary points and use the inpolygon function to 
%     reduce to the unique points within the boundary

[Xg,Yg]=meshgrid(min(MopXutmBound):max(MopXutmBound),...
    min(MopYutmBound):max(MopYutmBound));
% figure;plot(MopXutmBound,MopYutmBound,'k-');
% hold on;
% plot(Xg(:),Yg(:),'r.')
idx=find(inpolygon(Xg(:),Yg(:),MopXutmBound,MopYutmBound));
% plot(Xg(idx),Yg(idx),'c.')
% redefine unque x,y to based on inpolygon results
JumboXunq=Xg(idx);
JumboYunq=Yg(idx);

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
         [double(JumboYunq),double(JumboXunq)],...
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
        % Transect.Xutm=JumboXunq(idx);
        % Transect.Yutm=JumboYunq(idx);
        % Transect.UniqueXYindexes=idx;
        % hold on;plot(JumboXunq(idx),JumboYunq(idx),'.');
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
         MopTransectArea(n).Xutms=JumboXunq(idx);
         MopTransectArea(n).Yutms=JumboYunq(idx);
         MopTransectArea(n).UniqueXYindexes=idx;
         % hold on;plot(JumboXunq(idx),JumboYunq(idx),'.');
         % plot(MopTransectArea(n).Xutm,MopTransectArea(n).Yutm,'k.')
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
VolsDeep=[];
VolvDeep=[];
VolTot=[];
VolDeep=[];
VolShal=[];
VolvShal=[];
ns=0;
%figure
for s=ldx
    
 % check if profile has at least on epoint above MSL and one below -6m   
    igood=find(~isnan(Zf(s,:)));
    if max(Zf(s,:)) > VolumeElevCutoff && min(Zf(s,:)) < -6
        %clf
    %plot(X1Dmop(igood),Zf(s,igood));hold on;plot(X1Dmop,ZminGlobal)
    %plot(X1Dmop,X1Dmop*0+VolumeElevCutoff,'k:')
    
    if igood(end) < numel(X1Dmop)
    z1d=interp1([X1Dmop(igood) X1Dmop(end)],[Zf(s,igood) ZminGlobal(end)],...
        X1Dmop(igood(1)):X1Dmop(end),'linear','extrap');
    else
    z1d=interp1(X1Dmop(igood),Zf(s,igood),X1Dmop(igood(1)):X1Dmop(end),'linear','extrap');
    end
    xgood=X1Dmop(igood(1)):X1Dmop(end);
    z1dmin=ZminGlobal(igood(1):numel(X1Dmop));
    a1dgood=A1D(igood(1):numel(X1Dmop));
    %plot(xgood,z1d,'--')
    %pause
    
    ivol=find(z1d <= VolumeElevCutoff);%1.344);
    ivshal=find(z1d <= VolumeElevCutoff & z1dmin >= -2);
    ishal=find(z1d < -2 & z1dmin >= -4);
    ideep=find(z1dmin < -4 & z1dmin >= -6);
    ivdeep=find(z1dmin < -6 & z1dmin >= -8);
    isdeep=find(z1dmin < -8 & z1dmin >= -10);
    itot=find(z1d <= VolumeElevCutoff & z1dmin >= -10);
    % if ~isempty(ivol)
          ns=ns+1;
          Sdate(ns)=datetime(SA(s).Datenum,'ConvertFrom','datenum');
    %     VolCutElev(ns)=sum((z1d(ivol)-VolumeElevCutoff).*A1D(ia1d(ivol)));
          VolTot(ns)=sum((z1d(itot)-z1dmin(itot)).*a1dgood(itot));
          VolsDeep(ns)=sum((z1d(isdeep)-z1dmin(isdeep)).*a1dgood(isdeep));
          VolvDeep(ns)=sum((z1d(ivdeep)-z1dmin(ivdeep)).*a1dgood(ivdeep));
          VolDeep(ns)=sum((z1d(ideep)-z1dmin(ideep)).*a1dgood(ideep));
          VolShal(ns)=sum((z1d(ishal)-z1dmin(ishal)).*a1dgood(ishal));
          VolvShal(ns)=sum((z1d(ivshal)-z1dmin(ivshal)).*a1dgood(ivshal));
    % end
    end
end

Sdatenum=datenum(Sdate);
%numel(find(~isnan(VolTot)))
fprintf('\n%i valid volume estimates.\n\n',numel(find(~isnan(VolTot))));

TransectVols(Nmop).Mop=CurrentMopNumber;
TransectVols(Nmop).Vdatetimes=datetime([],[],[]);
TransectVols(Nmop).Vdatenums=[];
TransectVols(Nmop).VolTot=[];
TransectVols(Nmop).VolvShal=[];
TransectVols(Nmop).VolShal=[];
TransectVols(Nmop).VolDeep=[];
TransectVols(Nmop).VolvDeep=[];
TransectVols(Nmop).VolsDeep=[];

if numel(Sdatenum) > 0

TransectVols(Nmop).Vdatetimes=Sdate;
TransectVols(Nmop).Vdatenums=Sdatenum;
TransectVols(Nmop).VolTot=VolTot;
TransectVols(Nmop).VolvShal=VolvShal;
TransectVols(Nmop).VolShal=VolShal;
TransectVols(Nmop).VolDeep=VolDeep;
TransectVols(Nmop).VolvDeep=VolvDeep;
TransectVols(Nmop).VolsDeep=VolsDeep;

% if numel([TransectVols(Nmop).Vdatenums]) < numel([TransectVols(Nmop).VolTot])
%     numel([TransectVols(Nmop).Vdatenums]) 
%     numel([TransectVols(Nmop).VolTot])
%     fprintf('pausing..\n')
%     pause
% end

end

end

save MopSubaqueousVolumesEduardo.mat TransectVols 

% idx=find(~isoutlier(VolMin));
% figure('position',[100 100 1000 600]);
% plot(Sdate(idx),VolMin(idx)/100,'k.-');hold on;
% plot(Sdate(idx),VolDeep(idx)/100,'b.-')
% plot(Sdate(idx),VolShal(idx)/100,'r.-')
% plot(Sdate(idx),VolMin(idx)/100,'g.-')
% fprintf('%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n',mean(VolMin(idx)/100),...
% std(VolMin(idx)/100),...
% mean(VolShal(idx)/100),...
% std(VolShal(idx)/100),...
% mean(VolDeep(idx)/100),...
% std(VolDeep(idx)/100),...
% mean(VolMin(idx)/100),...
% std(VolMin(idx)/100));
% 
% %% --------------------------------------------------------------------
% %% Use accumarray to get area statistics
% 
% % JumboXutm=vertcat(SA(ldx).X);
% % JumboYutm=vertcat(SA(ldx).Y);
% % JumboZnavd88=vertcat(SA(ldx).Z);
% % reduce to MSL zone
% % imsl=find(JumboZnavd88 > 0.75 & JumboZnavd88 < 1.0);
% % JumboXutm=JumboXutm(imsl);
% % JumboYutm=JumboYutm(imsl);
% % JumboZnavd88=JumboZnavd88(imsl);
% 
% % bin all survey data points by placing in unique
% %  x,y data array
% [ux, ~, xidx] = unique(JumboXutm);
% [uy, ~, yidx] = unique(JumboYutm);
% 
% % array of counts of the number of points at each unique x/y combination
% zcount = accumarray([xidx(:), yidx(:)], 1);  
% %array of average z that fall into each unique x/y combination
% zavg = accumarray([xidx(:), yidx(:)], JumboZnavd88.')./zcount;
% zmedian = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @median); % 
% zmode = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @mode); % 
% zstd = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @std); % 
% zmin = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @min); % 
% zmax = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @max); % 
% zq75 = accumarray([xidx(:), yidx(:)], JumboZnavd88.',[], @(x) quantile(x,.75));
% %create a list of the z that fall into each unique x/y combination
% %zs = accumarray([xidx(:), yidx(:)],JumboZnavd88.', [], @(V) {V}, {});
% 
% % reduce arrays to 1d vectors of x,y points with accumarray z data 
% ii=isnan(zavg(:)) == 0; % 1d indices of valid data
% [i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% % final data vectors
% xutm=ux(i);yutm=uy(j); % utm points with data
% zavg=zavg(ii);
% zcount=zcount(ii);
% zmedian=zmedian(ii);
% zmin=zmin(ii);
% zmax=zmax(ii);
% zstd=zstd(ii);
% zmode=zmode(ii);
% zq75=zq75(ii);
% %%
% 
% vax2=axes('position',[.05 .05 .25 .4]);
% fill(MopLonBoundary,MopLatBoundary,[.95 .95 .95],'FaceAlpha',.2);
% hold on;
% plot([Mop.BackLon(MopNumber) Mop.OffLon(MopNumber)],...
%     [Mop.BackLat(MopNumber) Mop.OffLat(MopNumber)],'m*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber+1) Mop.OffLon(MopNumber+1)],...
%     [Mop.BackLat(MopNumber+1) Mop.OffLat(MopNumber+1)],'y*-','linewidth',2);
% 
% plot([Mop.BackLon(MopNumber-1) Mop.OffLon(MopNumber-1)],...
%     [Mop.BackLat(MopNumber-1) Mop.OffLat(MopNumber-1)],'y*-','linewidth',2);
% 
% plot_google_map('MapType', 'satellite','Alpha', 1);
% set(gca,'FontSize',12)
% 
% set(gca,'xlim',[Mop.BackLon(MopNumber)-0.0015 Mop.BackLon(MopNumber)+0.0015]);
% set(gca,'ylim',[Mop.BackLat(MopNumber)-0.0015 Mop.BackLat(MopNumber)+0.0015])
% 
% 
% % convert to lat lon
% [lat,lon]=utm2deg(xutm,yutm,repmat('11 S',[length(xutm) 1]));
% 
% % 
% %% plot min surface
% z=zmin;
% % plot colored points
% minz=-0.5;%min(z);%quantile(z,.05);
% maxz=5;%max(z);% zmax=quantile(z,.95);
% zrange=maxz-minz; % set max slope for coloring
% cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
% zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
% idx=find(~isnan(zscaled)); % non NaN points 
% scp=scatter(lon(idx), lat(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
% cb=colorbar;cb.Label.String=' Elevation (m, NAVD88) ';
% set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
% set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
% set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
% set(gca,'xtick',[],'ytick',[])
% plot_google_map('MapType', 'satellite','Alpha', 1);
% 
% title('Mobile LiDAR Global Minimum Elevation Surface','fontsize',14)
% 
% %% ------------------------------------------------------------------
% %%  Loop through individual surveys and interpolate-extrapolate in
% %%  2D to fill the subaerial mop area. 
% 
% 
% ns=0;
% Sdate=datetime([],[],[]);
% VolCutElev=[];
% VolMin=[];
% VolN=[];
% VolTot=[];
% 
% for s=ldx
%     % find indices of data points for this survey
%     %idx=find(SA(s).Datenum == JumboDatenum);
% 
%     % if more than one lidar survey on same day, use the last one
%     idx=find(s == JumboSurvnum);
%     %fprintf('%i\n',s)
%     %numel(idx)
%     % create scattered interpolant function with survey data
%     if numel(idx) > 3
%     F1 = scatteredInterpolant(JumboXutm(idx),JumboYutm(idx),...
%         JumboZnavd88(idx),'linear','linear');
%     % interpolate-extrapolate to get elevations at all points within
%     %  fixed subaerial mop area boundary
%     Zarea = F1(JumboXunq,JumboYunq);
%     else
%         Zarea=[];
%     end
% %% ------------------------------------------------------------------
% %% calculate the mop subaerial area volume above specified contour elevation
% %%
% 
% if numel(Zarea) > 0
%     ivol=find(Zarea >= VolumeElevCutoff);%1.344);
%     itot=find(Zarea > zmin);
%     if ~isempty(ivol)
%         ns=ns+1;
%         Sdate(ns)=datetime(SA(s).Datenum,'ConvertFrom','datenum');
%         VolCutElev(ns)=sum(Zarea(ivol)-VolumeElevCutoff);
%         VolMin(ns)=sum(Zarea(ivol)-zmin(ivol));
%         VolN(ns)=numel(find(JumboZnavd88(idx) >= VolumeElevCutoff))/numel(ivol);%numel(idx);
%         VolTot(ns)=sum(Zarea(itot)-zmin(itot));
%     end
% end
% 
% end
% % 
% % remove low lidar coverage cases
% % VolCutElev(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;
% % VolMin(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;
% % VolTot(VolN < MinFractionalCoverageAboveElevCutoff)=NaN;
% % 
% % 
% % vax3=axes('position',[.4 .55 .5 .35]);
% % plot using datenum so ginput function works on it
% % Sdatenum=datenum(Sdate);
% % numel(find(~isnan(VolTot)))
% % fprintf('\n%i valid volume estimates.\n\n',numel(find(~isnan(VolTot))));
% % 
% % LidarVols(Nmop).Mop=CurrentMopNumber;
% % LidarVols(Nmop).CutElev=VolumeElevCutoff;
% % LidarVols(Nmop).TotAreaMinElev=[];
% % LidarVols(Nmop).Vdatetimes=datetime([],[],[]);
% % LidarVols(Nmop).Vdatenums=[];
% % LidarVols(Nmop).VolCutElev=[];
% % LidarVols(Nmop).VolMin=[];
% % LidarVols(Nmop).VolTot=[];
% % 
% % 
% % if numel(Sdatenum) > 0
% % 
% % LidarVols(Nmop).Vdatetimes=Sdate;
% % LidarVols(Nmop).TotAreaMinElev=min(zmin(:));
% % LidarVols(Nmop).Vdatetimes=Sdate;
% % LidarVols(Nmop).Vdatenums=Sdatenum;
% % LidarVols(Nmop).VolCutElev=VolCutElev;
% % LidarVols(Nmop).VolMin=VolMin;
% % LidarVols(Nmop).VolTot=VolTot;
% % 
% % end
% 
% end
% 
% %save MopLidarVolumesMSLib.mat LidarVols 
% 
% %%
% %JumboVols=vertcat(LidarVols.VolMin);
% 
% 
% % p3=plot(Sdatenum,VolTot/1000,'g*-','linewidth',2);hold on
% % p2=plot(Sdatenum,VolMin/1000,'r*-','linewidth',2);hold on
% % p1=plot(Sdatenum,VolCutElev/1000,'k*-','linewidth',2);hold on
% % 
% % set(gca,'xtick',datenum(year(Sdate(1)):(year(Sdate(end))+1),1,1))
% % datetick('x','yyyy','keepticks')
% % set(gca,'fontsize',14);grid on;
% % ylabel('MOP Area Volume (K m^{3})')
% % legend([p1 p2 p3],['Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88'],...
% %     ['Vol above Minimum Elev Surface, Landward of ' num2str(VolumeElevCutoff,'%5.3fm') ' contour'],...
% %     'Total Vol above Minimum Surface','location','north')
% % yl=get(gca,'ylim');set(gca,'ylim',[yl(1) yl(1)+(yl(2)-yl(1))*1.2]);
% % title(['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'])
% % 
% % zoom on;
% % option to edit a survey
% % Vedit=uicontrol(gcf,'Style','Pushbutton','Position',[580 650 150 25],...
% % 'Callback','GetVolSurveyNumber;EditSA','String','Edit a Survey'); 
% % set(Vedit,'backgroundcolor','yellow')
% 
% 
% %%
% %vax4=axes('position',[.4 .08 .55 .35]);
% % sy=year(Sdate);sm=month(Sdate);
% % hold on;
% % for y=min(sy):max(sy)
% %     if y == max(sy)
% %      plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'m*-',...
% %      'linewidth',3,'displayname',num2str(y));hold on
% %     else       
% %      plot(month(Sdate(sy == y))+day(Sdate(sy == y))/31,VolCutElev(sy == y)/1000,'*-',...
% %      'linewidth',2,'displayname',num2str(y));hold on
% %     end
% % end
% % set(gca,'xlim',[1 13],'xtick',1:12,'fontsize',14);grid on;legend('location','eastoutside')
% % xlabel('Month of Year');ylabel('MOP Area Volume (K m^{3})')
% % title({['MOP ' num2str(MopNumber) ' Mobile LIDAR Subaerial Volumes'],...
% %     [' Seasonal view, Vol above ' num2str(VolumeElevCutoff,'%5.3fm') ' NAVD88']})
% 
