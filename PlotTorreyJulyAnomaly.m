% north torrey july sand elevation anomaly figure

% torrey
MopStart=find(strcmp([Mop.Name],'D0580'));
MopEnd=find(strcmp([Mop.Name],'D0594'));

DateStart=datenum(2024,7,6);
DateEnd=datenum(2024,7,9);

%CS=SAcombineMops(MopStart,MopEnd);
%CS=SGcombineMops(MopStart,MopEnd);
%% build Global month mean grid for mop range
Mon=7;
CGM=GMcombineMopsMonth(MopStart,MopEnd,Mon);

xmin=min(vertcat(CGM.X2D));xmax=max(vertcat(CGM.X2D));
ymin=min(vertcat(CGM.Y2D));ymax=max(vertcat(CGM.Y2D));
xmin=floor(xmin);xmax=ceil(xmax);
ymin=floor(ymin);ymax=ceil(ymax);
x=CGM.X2D;
y=CGM.Y2D;
% x,y grid arrays
[Xg,Yg]=meshgrid(xmin:xmax,ymin:ymax);
idx=sub2ind(size(Xg),y-ymin+1,x-xmin+1);
% Global Month Mean z grid array
Zmean=Xg*NaN; % initialize as NaNs
Zmean(idx)=CGM.Z2Dmean; % assign valid z grid data points

%% Build recent July grid  

CS=SGcombineMops(MopStart,MopEnd);

Zjuly=Xg*NaN; % initialize as NaNs, same grid as global mean
%x=[];y=[];z=[];
% first get gridded jumbo
CS=SAcombineMops(MopStart,MopEnd);
jdx=find( strcmp({CS.Source},'Gps') &...
    [CS.Datenum] >= DateStart & [CS.Datenum] <= DateEnd);
% then beach lidar
ndx=find((strcmp({CS.Source},'AtvMR') | strcmp({CS.Source},'Trk')) &...
    [CS.Datenum] >= DateStart & [CS.Datenum] <= DateEnd);

for idx = [jdx ndx] 
x=CS(idx).X;
y=CS(idx).Y;
z=CS(idx).Z;
gdx=sub2ind(size(Xg),y-ymin+1,x-xmin+1);
% overlay gridded recent survey data on grid
Zjuly(gdx)=z; % assign valid z grid data points
end

Anom=Zjuly-Zmean;
mdx=find(~isnan(Anom(:)));
%figure;surf(Xg,Yg,Zjuly-Zmean);colormap(flipud(polarmap));shading flat;view(2)
%%
figure('position',[1          12        400         783]);
[lat,lon]=utm2deg(Xg(mdx),Yg(mdx),repmat('11 S',[length(Xg(mdx)) 1]));
[ScatterPlot,ColorBarPlot]=ColorScatterPolarmap(lon,lat,Anom(mdx));
hold on
load MopTableUTM.mat
% for n=MopStart:MopEnd
% plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
% text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
% end
%plot(lon2(idx2),lat2(idx2),'m.')
plot_google_map('MapType', 'satellite')
% set(gca,'clim',[-13 -4])
% set(gca,'clim',[-1 6])
set(gca,'fontsize',16)

% title({[datestr(CS(idx).Datenum) ' |  Multibeam '],...
%     'Mops 645 to 655'},...
%     'fontsize',16);
title({'Penasquitos Inlet, Torrey Pines','July 2024 Sand Elev. Anomaly'},'fontsize',18)
set(gca,'xtick',[],'ytick',[]);box on;
pos=get(gca,'position');
set(gca,'position',[pos(1)-0.085 0.05 pos(3)-.06 0.85])
set(gca,'xlim',[-117.2684 -117.2586],'ylim',[32.9240   32.9400]);

makepng('TorreyJulyElevationAnomaly2024.png')
%ScatterPlotBeachUTM([CGM.X2D],[CGM.Y2D],[CGM.Z2Dmean],'2d');
% combine load global mean files