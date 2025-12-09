% Example code to plot a mutibeam survey on a google map

% takes a few minutes to run

%% settings -----------------------------------------------------------

% path to matlab functions 
addpath '/Volumes/group/MOPS'
addpath '/Volumes/group/MOPS/toolbox'
% full filename of multibeam las file to read
dataFilename='/Volumes/group/Multibeam/20241023_07278_07288_Seacliff_multibeam/20241022-23_Seacliff_preliminary.las';
% Also need to know the utm zone of the data for conversion to lat lon
% Using zone 11 S for Seacliff to get it to plot correctly even though
% it is technically in zone 10 S.
utmzone='11 S';
% and finally the Mop range of transect lines to plot with data
MopStart=7278;
MopEnd=7288;

%% -----------------------------------------------------------------------

load MopTableUTM.mat

%% 1. read in the .las multibeam survey to struct array d
d=lasdata(dataFilename,'loadall');

%% 2. extract xutm, yutm and z of survey points
    xutm=d.x;
    yutm=d.y;-*
    z=d.z;

%% 3. reduce to 1m utm spatial average values
%----------------------------------------
% reduce to 1m spatial averages  
%----------------------------------------
 
Res=1; % 1m spatial resolution

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

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final shore box data vectors
xutm=ux(i);yutm=uy(j);
zavg=zavg(ii);

fprintf(1,...
    'Reduced to %g , %g x %g meter spatially averaged survey points\n',...
    length(xutm),Res,Res);

%% 4. transform to lat lons based on UTM zone
[lat,lon]=utm2deg(xutm,yutm,repmat(utmzone,[length(xutm) 1]));

%% 5. plot lat,lon, depths as colored dots on google map
figure('position',[1          12        1275         783]);%plot(lon,lat,'y.');
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
hold on;
[ScatterPlot,ColorBarPlot]=ColorScatter2d(lon,lat,zavg);
hold on
load MopTableUTM.mat
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
plot_google_map('MapType', 'satellite')
set(gca,'fontsize',16)
title('Seacliff State Beach Shallow Multibeam Survey')
xlabel('Longitude');ylabel('Latitude');

% make a png file of the figure
pfile='SeacliffMultibeam.png';
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','-r300','-loose',pfile);
fprintf(1,'Output image written to: %s\n',pfile)
