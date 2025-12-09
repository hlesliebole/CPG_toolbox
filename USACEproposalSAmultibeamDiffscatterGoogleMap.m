clearvars
load MopTableUTM.mat

% MopStart=find(strcmp([Mop.Name],'D0645'));
% MopEnd=find(strcmp([Mop.Name],'D0654'));

MopStart=find(strcmp([Mop.Name],'D0649'));
MopEnd=find(strcmp([Mop.Name],'D0657'));

CS=SAcombineMops(MopStart,MopEnd);
idx=find(strcmp({CS.Source},'Multibeam'));
%idx=find([CS.Datenum] == datenum(2023,11,21));
%idx=idx(2);
%idx=idx(2);

x1=vertcat(CS(idx(1)).X);
y1=vertcat(CS(idx(1)).Y);
z1=vertcat(CS(idx(1)).Z);

x2=vertcat(CS(idx(3)).X);
y2=vertcat(CS(idx(3)).Y);
z2=vertcat(CS(idx(3)).Z);

% find common grid points
[l,li]=ismember([x1 y1],[x2 y2],'rows');
% depth differences
zdiff=z2(li(l > 0))-z1(l > 0);
bias=mean(zdiff);
zdiff=zdiff-bias;
% 
% % trim higher elevations
% zidx=find(z < 6);
% x=x(zidx);y=y(zidx);z=z(zidx);
% convert to lat lon
%[lat,lon]=utm2deg(x1(l >0),y1(l > 0),repmat('11 S',[length(x1(l > 0)) 1]));
[lat,lon]=utm2deg(x2(li(l >0)),y2(li(l > 0)),repmat('11 S',[length(x2(li(l > 0))) 1]));

%figure;
figure('position',[1    12   867   783]);%plot(lon,lat,'y.');
% for n=MopStart:MopEnd
% plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
% text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
% end
hold on;
[ScatterPlot,ColorBarPlot]=ColorScatter2d(lon,lat,zdiff);
hold on
% load MopTableUTM.mat
% for n=MopStart:MopEnd
% plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
% text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
% end
%plot(lon2(idx2),lat2(idx2),'m.')
set(gca,'xlim',[-117.2851 -117.2738]);
set(gca,'ylim',[32.9847   32.9935]);
%%
%delete(pl)
%delete(tl)
pl=plot([-117.281 -117.281+cosd(33)*.5*(0.01/1.11)],[32.985 32.985],'w-','linewidth',3);
tl=text(-117.2796,32.9852,'500 m','fontsize',22,'fontweight','bold','color','w');
%%
% delete(pd)
% delete(td)
% delete(td2)
delete(ps)
xd=[-117.2778 -117.2750];
yd=[32.9907 32.9912];
pd=plot(xd,yd,'w--','linewidth',3);
td=text(-117.2774,32.9904,'DREDGE PIPE','color','w','fontweight','bold',...
    'fontsize',20,'rotation',13);
td2=text(-117.2776,32.9918,{'FLETCHER','    COVE'},'color','w','fontweight','bold',...
    'fontsize',20,'rotation',0);
xs =[-117.2811 -117.2778];
ys =[32.9917  32.9895];
ps=plot([xs(1) xs(2) xs(2) xs(1) xs(1)],...
    [ys(1) ys(1) ys(2) ys(2) ys(1)],'-','Color',[.6 .6 .6],'linewidth',8);


%%

plot_google_map('MapType', 'satellite')
%set(gca,'clim',[-4 4])
set(gca,'fontsize',16)
ColorBarPlot.FontSize=22;
ColorBarPlot.Label.String='Elevation Change (m)';
title({'Solana Beach: Elevation change 21 Nov 2023',...
    '(pre-nourishment) to 1 Mar 2024 (during nourishment)'},...
    'fontsize',24)
set(gca,'xtick',[],'ytick',[])
set(gcf,'InvertHardcopy','off')
makejpg('MultibeamDredgePipe.jpg')
% SA=CS;
% jumbo=find(contains({SA.File},'umbo') & [SA.Datenum] == datenum(2023,10,17));
% SA=SA(jumbo(end));
% 
% %% find intersection of 2 data sets
% % idx=find(ismember([SA.X SA.Y],[x2 y2],'rows'));
% % ndx=find(ismember([x2 y2],[SA.X SA.Y],'rows'));
% 
% %% plot differneces on map
% 
% %% make scatterplot of depths
% %zavg=zsave;
% %bdx=find(abs(SA.Z(idx)-z2(ndx)) > 1);z2(ndx(bdx))=NaN;
% z2=z2(li(l > 0));z1=z1(l > 0);
% 
% rmse=sqrt(mean((z2-z1).^2,'omitnan'));
% bias=mean(z1-z2);
% brrmse=sqrt(mean(( (z1-bias)-z2).^2,'omitnan'));
% N=numel(z1);
% figure('position',[380   178   707   564]);
% plot(z2,z1,'k.');hold on;plot([-15 0],[-15 0],'m-','linewidth',2)
% grid on;set(gca,'xlim',[-15 0],'ylim',[-15 0],'fontsize',14);
% xlabel([datestr(CS(idx(3)).Datenum) ' Multibeam (m, NAVD88)']);
% ylabel([datestr(CS(idx(1)).Datenum) ' Multibeam Depth (m, NAVD88)']);
% title({'Solana Multibean vs Multibeam :',['Comparison of N= ' num2str(N) ' overlapping 1m spatially-averaged depths']},'fontsize',16)
% text(-14,-2,['RMSE = ' num2str(rmse,'%6.3fm')],'fontsize',16)
% text(-14,-3,['BIAS = ' num2str(bias,'%6.3fm')],'fontsize',16)
% text(-14,-4,['RMSE (BIAS removed) = ' num2str(brrmse,'%6.3fm')],'fontsize',16)
% set(gcf,'inverthardcopy','off')
%makepng('MultibeamVsMultibeamSolanaJan2024.png')


function [ScatterPlot,ColorBarPlot]=ColorScatter2d(x,y,z)

%zmin=quantile(z,.05);zmax=quantile(z,.95);
zmin=-0.5;zmax=0.5;
%zmin=-15;zmax=-3;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
cm =flipud(polarmap(64));
%cm =jet(64);

scp=scatter(x(idx), y(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
%view(-10,80)
colormap(flipud(polarmap(64)))
%colormap(jet(64))

cb=colorbar;cb.Label.String='Elevation Change (m, navd88)';
set(gca,'clim',[zmin zmax]);
% set(gca,'xlim',[min(x) max(x)]);
% set(gca,'ylim',[min(y) max(y)]);
%set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end