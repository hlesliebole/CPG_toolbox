% code to build Solana prenoursihment baseline grid

clearvars
close all

load SolanaShoreboxMap.mat
figure('position',[ 33          53        1381         739]);
image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');hold on;
set(gca,'dataaspectratio',[1 1 1],'xdir','normal','fontsize',16);
set(gca,'ylim',[min(ygimg) max(ygimg)],'xlim',[min(xgimg) max(xgimg)])
pos=get(gca,'position');
%set(gca,'position',[.5 pos(2) pos(3) pos(4)],'fontsize',14,'linewidth',2)
xlabel('Alongshore Distance (m)','fontsize',18);
ylabel('Cross-shore Distance (m)','fontsize',18);
title('Solana Beach Nourishment Pad Evolution :  10 May vs 16 Apr 2024','fontsize',22)
xl=get(gca,'xlim');
hold on;
for n=3:31
%plot([MopSB.BackLon(n) MopSB.OffLon(n)],[MopSB.BackLat(n) MopSB.OffLat(n)],'m-','linewidth',2)
%text(MopSB.OffLat(n)+300,MopSB.OffLon(n),[ '   ' MopSB.Name{n}],'color','y','fontsize',16,'fontweight','bold')
end
for n=4:2:30
%plot([MopSB.BackLat(n) MopSB.OffLat(n)],[MopSB.BackLon(n) MopSB.OffLon(n)],'m-','linewidth',2)
% plot(MopSB.OffLon(n),MopSB.OffLat(n),'y.','markersize',15)
% text(MopSB.OffLon(n),MopSB.OffLat(n),[ '  ' MopSB.Name{n}],'color','y','fontsize',18,'fontweight','bold','rotation',90)
end

tf=text(1850,-100,{'Fletcher','  Cove'},'color','w','fontsize',22,'fontweight','bold');
set(gcf,'inverthardcopy','off','color','w');
%makepng('SolanaShoreBox.png')

%% Build a shorebox baseline elev grid 
% make 2d X,Y grid arrays
[X,Y]=meshgrid(xgimg,ygimg);

load SolanaBaselineGrid
load SolanaPrenourishmentGrid.mat
Z0P=Z0;clear Z0;
load SolanaPostnourishmentGrid.mat  % 16 Apr 24
Z0=Z;clear Z;
load SolanaPostnourishmentGrid26Apr24.mat
Z1=Z;clear Z;
load SolanaPostnourishmentGrid01May24.mat
Z2=Z;clear Z;
load SolanaPostnourishmentGrid10May24.mat
Z3=Z;clear Z;
load SolanaPostnourishmentGrid24May24.mat
Z3=Z;clear Z;

load SolanaShoreboxMap.mat
idx=find(Z0(:) > -12.5 & Z0(:) < 6 & ~isnan(Z0(:)) & ~isnan(Z1(:)) & ~isnan(Z2(:)) & ~isnan(Z3(:)));
dZ1=Z0*NaN;dZ1(idx)=Z1(idx)-Z0(idx);
dZ2=Z0*NaN;dZ2(idx)=Z2(idx)-Z0(idx);
dZ3=Z0*NaN;dZ3(idx)=Z3(idx)-Z0(idx);

hold on;
imagesc(xgimg,ygimg,dZ3,'AlphaData',~isnan(dZ3));hold on;
set(gca,'clim',[-1 1]);colormap(flipud(polarmap));cb=colorbar;
cb.Label.String='Elevation Change (m)';
 cb.Ticks=[-1:.2:1];
set(gca,'ylim',[-200 850]);
set(gca,'position',[0.1259    0.3100    0.7503    0.8150]);

% [C,h] =contour(xgimg,ygimg,Z0P-0.775,[0 -4 -10],'m:','linewidth',2);
% text(-275,500,{'PreNourishment',' 0, -4, -10m MSL'},'color','m','fontsize',18,'fontweight','bold')
set(gca,'ylim',[-100 250]);
set(gca,'dataaspectratio',[1 .4 1]);hold on;
plot([400 400],[-100 100],'y-','linewidth',5);plot([800 800],[-100 100],'y-','linewidth',5);
plot([2700 2700],[-100 0],'y-','linewidth',5);plot([3000 3000],[-100 0],'y-','linewidth',5);
text(600,-50,{'South of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(2850,-50,{'North of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(1750,-30,{'Pad Reach'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
%clabel(C,h)
% title({'Solana Beach Nourishment',[datestr(sdate) ' Truck LiDAR Survey'],...
%     'Observed Vol Change above MSL^{**} = 440K m^3 (525K planned) ** new subaqueous toe not included'})
%%
ax2=axes('position',[0.1259    0.3100    0.7503    0.1650]);

p1=plot(xgimg,sum(dZ1,'omitnan'),'-','linewidth',2);hold on;
p2=plot(xgimg,sum(dZ2,'omitnan'),'-','linewidth',2);hold on;
p2=plot(xgimg,sum(dZ3,'omitnan'),'-','linewidth',2);hold on;
plot(xgimg,0*sum(dZ3,'omitnan'),'k:','linewidth',1);
set(gca,'fontsize',16);grid on;set(gca,'ylim',[-30 30],'ytick',-20:10:30,'xlim',[xgimg(1) xgimg(end)])
xlabel(' ','fontsize',18);
ylabel({'Vol Change','(m^{3}/m-shore)'})

plot([400 400],[-30 30],'k:','linewidth',2);
plot([800 800],[-30 30],'k:','linewidth',2);
plot([2700 2700],[-30 30],'k:','linewidth',2);
plot([3000 3000],[-30 30],'k:','linewidth',2);

idx=find(xgimg >= 400 & xgimg < 800);
SP1=sprintf('%+8i',round(sum(sum(dZ1(:,idx),'omitnan'))));
SP2=sprintf('%+8i',round(sum(sum(dZ2(:,idx),'omitnan'))));
dSP2=sprintf('(%+8i)',round(sum(sum(dZ2(:,idx),'omitnan')))-round(sum(sum(dZ1(:,idx),'omitnan'))));
SP3=sprintf('%+8i',round(sum(sum(dZ3(:,idx),'omitnan'))));
dSP3=sprintf('(%+8i)',round(sum(sum(dZ3(:,idx),'omitnan')))-round(sum(sum(dZ2(:,idx),'omitnan'))));
% exclude > pad elevations from pad reach change calc
dZ1p=dZ1;dZ1p(Z1 > 4.4)=NaN;
dZ2p=dZ2;dZ2p(Z2 > 4.4)=NaN;
%%
nc=0;
idx=find(xgimg >= 800 & xgimg < 2700);
for zcut=1:.01:5
dZ1p=dZ3;dZ1p(Z3 > zcut)=NaN;
nc=nc+1;
zc(nc)=zcut;
vcut(nc)=round(sum(sum(dZ1p(:,idx),'omitnan')));
end

figure('position',[ 265   109   809   643]);
plot(zc,vcut,'b.-','linewidth',2);
hold on;grid on;
plot([4.3 4.3],[min(vcut) max(vcut)],'k--','linewidth',2);
p1=plot([4.7 4.7],[min(vcut) max(vcut)],'k--','linewidth',2);
set(gca,'fontsize',16);legend(p1,'Primary Terrace Elv Range','location','northeast')
xlabel('Pad Volume Change Estimate Elevation Cutoff (m, NAVD88)');
ylabel('Estimated 16 Apr to 24 May 2024 Pad (Beach Face) Volume Change');
title('Pad Volume Change Estimate Sensitivity to the Pad Terrace Elevation Cutoff');
makepng('SolanaPadElevationVolumeSensitivity24may.png');

%%
