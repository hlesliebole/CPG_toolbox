load MopTableUTM.mat

MopStart=find(strcmp([Mop.Name],'D0645'));
MopEnd=find(strcmp([Mop.Name],'D0657'));

MopStart=find(strcmp([Mop.Name],'D0632'));
MopEnd=find(strcmp([Mop.Name],'D0647'));

MopStart=find(strcmp([Mop.Name],'D0582'));
MopEnd=find(strcmp([Mop.Name],'D0590'));

CS=SAcombineMops(MopStart,MopEnd);
idx=find(strcmp({CS.Source},'Multibeam'));
%idx=find([CS.Datenum] == datenum(2023,11,21));
%idx=idx(2);
%idx=idx(end);
idx=idx(1);

x=vertcat(CS(idx).X);
y=vertcat(CS(idx).Y);
z=vertcat(CS(idx).Z);
% trim higher elevations
zidx=find(z < 6);
x=x(zidx);y=y(zidx);z=z(zidx);
% convert to lat lon
[lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));

%figure;
figure('position',[1          12        1275         783]);%plot(lon,lat,'y.');
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
hold on;
[ScatterPlot,ColorBarPlot]=ColorScatter2d(lon,lat,z);
%%
hold on
load MopTableUTM.mat
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
%plot(lon2(idx2),lat2(idx2),'m.')
plot_google_map('MapType', 'satellite')
%set(gca,'clim',[-13 -4])
set(gca,'clim',[-11 -1])
set(gca,'fontsize',16)

title({[datestr(CS(idx).Datenum) ' |  Multibeam '],...
    'Mops 582 to 590'},...
    'fontsize',16)
set(gcf,'inverthardcopy','off')
makepng(['Multibeam' datestr(CS(idx).Datenum,'YYYYmmDD') 'Mops' num2str(MopStart) 'to' num2str(MopEnd) '.png'])

mbDate=CS(idx).Datenum;
SA=CS;
x2=x;y2=y;z2=z;
jumbo=find(contains({SA.File},'umbo') & [SA.Datenum] == datenum(2024,1,9));
jumbo=find(contains({SA.File},'umbo') & [SA.Datenum] == datenum(2023,12,14));
jumbo=find(contains({SA.File},'umbo') & [SA.Datenum] == datenum(2025,1,28));
%jumbo=find(contains({SA.File},'umbo') & [SA.Datenum] == datenum(2023,10,17));
SA=SA(jumbo(end));

%% find intersection of 2 data sets
idx=find(ismember([SA.X SA.Y],[x2 y2],'rows'));
ndx=find(ismember([x2 y2],[SA.X SA.Y],'rows'));

%% plot differneces on map

%% make scatterplot of depths
%zavg=zsave;
bdx=find(abs(SA.Z(idx)-z2(ndx)) > 1);z2(ndx(bdx))=NaN;
rmse=sqrt(mean(( SA.Z(idx)-z2(ndx)).^2,'omitnan'));
bias=mean(z2(ndx)-SA.Z(idx),'omitnan');
brrmse=sqrt(mean(((z2(ndx)-bias)-SA.Z(idx)).^2,'omitnan'));
N=numel(idx)-numel(bdx);
figure('position',[380   178   707   564]);
plot(SA.Z(idx),z2(ndx),'k.');hold on;plot([-12 0],[-12 0],'k--')
grid on;set(gca,'xlim',[-12 0],'ylim',[-12 0],'fontsize',14);
xlabel([datestr(SA.Datenum) ' Jetski Depth (m, NAVD88)']);
ylabel([datestr(mbDate) ' Multibeam Depth (m, NAVD88)']);
%title({'Solana Multibeam vs Jeski :',['Comparison of N= ' num2str(N) ' overlapping 1m spatially-averaged depths']},'fontsize',16)
title({'Torrey Multibeam vs Jeski :',['Comparison of N= ' num2str(N) ' overlapping 1m spatially-averaged depths']},'fontsize',16)
text(-11,-2,['RMSE = ' num2str(rmse,'%6.3fm')],'fontsize',16)
text(-11,-3,['BIAS = ' num2str(bias,'%6.3fm')],'fontsize',16)
text(-11,-4,['RMSE (BIAS removed) = ' num2str(brrmse,'%6.3fm')],'fontsize',16)
%makepng('MultibeamVsJetskiSolanaJan2024.png')
%makepng('MultibeamVsJetskiSolanaNov2023.png')
%makepng('MultibeamVsJetskiSolana26Jan2024.png')
makepng('MultibeamVsJetskiTorrey31Jan2025.png')



function [ScatterPlot,ColorBarPlot]=ColorScatter2d(x,y,z)

%zmin=quantile(z,.05);zmax=quantile(z,.95);
zmin=-11;zmax=-1;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
cm =jet(64);

scp=scatter(x(idx), y(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
%view(-10,80)
colormap(jet(64))

cb=colorbar;cb.Label.String='Elevation (m, navd88)';
set(gca,'clim',[zmin zmax]);
% set(gca,'xlim',[min(x) max(x)]);
% set(gca,'ylim',[min(y) max(y)]);
%set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end