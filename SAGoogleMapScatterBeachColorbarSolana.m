load MopTableUTM.mat

MopStart=find(strcmp([Mop.Name],'D0645'));
MopEnd=find(strcmp([Mop.Name],'D0657'));
MopStart=find(strcmp([Mop.Name],'D0632'));
MopEnd=find(strcmp([Mop.Name],'D0647'));
MopStart=find(strcmp([Mop.Name],'D0636'));
MopEnd=find(strcmp([Mop.Name],'D0665'));

% san dieguito
MopStart=find(strcmp([Mop.Name],'D0634'));
MopEnd=find(strcmp([Mop.Name],'D0638'));

% solana nourishment
MopStart=find(strcmp([Mop.Name],'D0634'));
MopEnd=find(strcmp([Mop.Name],'D0667'));

DateStart=datenum(2024,7,6);
DateEnd=datenum(2024,7,9);

x=[];y=[];z=[];
% first get gridded jumbo
CS=SGcombineMops(MopStart,MopEnd);
ndx=find( strcmp({CS.Source},'Gps') &...
    [CS.Datenum] >= DateStart & [CS.Datenum] <= DateEnd);
for idx = ndx 
x=[x' vertcat(CS(idx).X)']';
y=[y' vertcat(CS(idx).Y)']';
z=[z' vertcat(CS(idx).Z)']';
end

% now overlay beach lidar
CS=SAcombineMops(MopStart,MopEnd);

ndx=find((strcmp({CS.Source},'AtvMR') | strcmp({CS.Source},'Trk')) &...
    [CS.Datenum] >= DateStart & [CS.Datenum] <= DateEnd);
%idx=find([CS.Datenum] == datenum(2023,11,21));
%idx=idx(2);
%idx=idx(end-4);
%idx=305;

for idx = ndx 
x=[x' vertcat(CS(idx).X)']';
y=[y' vertcat(CS(idx).Y)']';
z=[z' vertcat(CS(idx).Z)']';
end

% trim higher elevations
zidx=find(z < 10);
x=x(zidx);y=y(zidx);z=z(zidx);
% convert to lat lon
[lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));

%figure;

% san dieguito figure
%figure('position',[1          12        600         783]);%plot(lon,lat,'y.');

% solana figure
figure('position',[ 100    59   400   738]);

% for n=MopStart:MopEnd
% plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
% text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
% end
hold on;
%[ScatterPlot,ColorBarPlot]=ColorScatter2d(lon,lat,z);
%ScatterPlotBeachUTM(x,y,z,PlotMethod)
ScatterPlotBeachLatLon(lat,lon,z,'2d');
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

% title({'San Dieguito Inlet, Del Mar CA',datestr(CS(idx).Datenum)},'fontsize',18)
% set(gca,'xtick',[],'ytick',[]);box on;
% pos=get(gca,'position');
% set(gca,'position',[pos(1)-0.085 0.05 pos(3)-.12 0.85])

title({'Solana Beach, CA',datestr(CS(idx).Datenum)},'fontsize',18)
set(gca,'xtick',[],'ytick',[]);box on;
set(gca,'xlim',[ -117.2878 -117.2662],'ylim',[32.9700   33.0050]);
%%
pos=get(gca,'position');
set(gca,'position',[pos(1)-0.085 0.05 pos(3)-.12 0.85])
BeachColorbar
%%
makepng(['AtvMR' datestr(CS(idx).Datenum,'YYYYmmDD') 'Mops' num2str(MopStart) 'to' num2str(MopEnd) '.png'])