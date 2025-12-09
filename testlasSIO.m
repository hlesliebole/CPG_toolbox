
% %dfile='/volumes/drone/data/torrey/20220216/20220216_00581_00590_TorreyCobble_RTKdrone_epoch2010_geoid12b.las';
% dfile='/Users/William/Desktop/Scripps_Flight07_classified_NAVD88_corrected.las';
% d=lasdata(dfile,'loadall');
% x=d.x;
% y=d.y;
% z=d.z;
% delete(d);
% save SIOsurfzoneLiDAR.mat x y z



%% reduce to 1m avg data

load SIOsurfzoneLiDAR.mat x y z
xutm=x;yutm=y;
clear x y


% idx=find(z > -25.5);
% xm=median(xutm(idx));
% ym=median(yutm(idx));
% %idx=find(xutm > xm & xutm < xm+10 & yutm > ym & yutm < ym+10);
% xa=xutm(idx);ya=yutm(idx);za=z(idx);
% za(abs(diff(xa)) > 1)=NaN;za(abs(diff(ya)) > 1)=NaN;
% figure;
% plot3(xa,ya,za,'k-')
% %ScatterPlotSubaqueousUTMdetail(xa,ya,za,'3d')
% set(gca,'dataaspectratio',[1 1 1])

%----------------------------------------
% reduce to 1m spatial averages  
%----------------------------------------
 
Res=.1; % 1m spatial resolution

% round survey x,y to desired resolution
xr=Res*round(xutm/Res); % round to Res meter spatial resolution
yr=Res*round(yutm/Res); % 

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(xr);
[uy, ~, yidx] = unique(yr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average of z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
zmin = accumarray([xidx(:), yidx(:)], z.',[],@min);
zmax = accumarray([xidx(:), yidx(:)], z.',[],@max);
%cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
%tavg = accumarray([xidx(:), yidx(:)], t.')./zcount;
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final shore box data vectors
xutm=ux(i);yutm=uy(j);
zmin=zmin(ii);
zmax=zmax(ii);
% figure;ScatterPlotBeachUTM(xutm,yutm,zmax,'3d')
% hold on;
% ScatterPlotBeachUTM(xutm,yutm,zmin,'3d')
% set(gca,'dataaspectratio',[1 1 .2],'zlim',[-2 4])
% view(-37.5,5)
figure;plot3(xutm,yutm,zmax,'k.')
hold on;
plot3(xutm,yutm,zmin,'r.')
set(gca,'dataaspectratio',[1 1 .1],'zlim',[-2 4])
view(-37.5,5)
grid on;

%%
%cmode=cmode(ii);
%tavg=tavg(ii);

%fprintf('Min-Max Classification Mode: %d %d\n',min(cmode),max(cmode))

%--------------------------------------------------------
% option to remove spatially averaged data based on less 
%  than N survey points
%
% zcount=zcount(ii); 
% N=3;
%i=find(zcount < N);xutm(i)=[];yutm(i)=[];zavg(i)=[];
%  zcount(i)=[];cmode(i)=[];
%--------------------------------------------------------

%end % end if for npts > 0

fprintf(1,...
    'Reduced to %g , %g x %g meter spatially averaged survey points\n',...
    length(xutm),Res,Res);

% multibeam 1m data data is now in xutm yutm zavg

%% get combine SA data for mops 519-527
SA=SAcombineMops(509,512);
SA=SA(122);

%% find intersection of 2 data sets
idx=find(ismember([SA.X SA.Y],[xutm yutm],'rows'));
ndx=find(ismember([xutm yutm],[SA.X SA.Y],'rows'));

%% plot differneces on map

%% make scatterplot of depths
%zavg=zsave;
bdx=find(abs(SA.Z(idx)-zavg(ndx)) > 1);zavg(ndx(bdx))=NaN;
rmse=sqrt(mean(( SA.Z(idx)-zavg(ndx)).^2,'omitnan'));
N=numel(idx)-numel(bdx);
figure('position',[380   178   707   564]);
plot(SA.Z(idx),zavg(ndx),'k.');hold on;plot([-25 0],[-25 0],'k--')
grid on;set(gca,'xlim',[-25 0],'ylim',[-25 0],'fontsize',14);
xlabel('Jetski Depth (m, NAVD88)');
ylabel('LiDAR Depth (m, NAVD88)');
title(['Comparison of N= ' num2str(N) ' overlapping 1m spatially-averaged depths'],'fontsize',16)
text(-24,-3,['RMSE = ' num2str(rmse,'%6.3fm')],'fontsize',16)
makepng('SIOsurfzoneLiDAR.png')

%%
% 
% figure;
% ScatterPlotSubaqueousUTM(x(1:10:end),y(1:10:end),z(1:10:end),-25,'2d')
% 
% 
% idx=find(x > 476025 & x < 476323 & y > 3637346 & y < 3637581);
% xc=x(idx);yc=y(idx);zc=z(idx);
% figure;
% ScatterPlotSubaqueousUTM(xc,yc,zc,-25,'2d')
% 
% [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));
% 
% figure('position',[219          14        1023         783]);
% ScatterPlotSubaqueous(lon(1:10:end),lat(1:10:end),z(1:10:end),-25,'2d')
% ax1=gca;set(gca,'fontsize',14)
% plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1);
% for m=519:527
% PlotLabelMopTransect(m,'2d','w','ShadeOff');
% end
% title('SIO CPG Shallow Water Multibeam','fontsize',18)