% code to build Solana prenoursihment baseline grid

clearvars
close all

load SolanaShoreboxMap.mat
figure('position',[ 33          53        1381         739]);
ax1=axes;
image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');hold on;
set(gca,'dataaspectratio',[1 1 1],'xdir','normal','fontsize',16);
set(gca,'ylim',[min(ygimg) max(ygimg)],'xlim',[min(xgimg) max(xgimg)])
pos=get(gca,'position');
%set(gca,'position',[.5 pos(2) pos(3) pos(4)],'fontsize',14,'linewidth',2)
xlabel('Alongshore Distance (m)','fontsize',18);
ylabel('Cross-shore Distance (m)','fontsize',18);

xl=get(gca,'xlim');
hold on;
for n=3:31
plot([MopSB.BackLon(n) MopSB.OffLon(n)],[MopSB.BackLat(n) MopSB.OffLat(n)],'m-','linewidth',2)
%text(MopSB.OffLat(n)+300,MopSB.OffLon(n),[ '   ' MopSB.Name{n}],'color','y','fontsize',16,'fontweight','bold')
end
for n=4:2:30
% plot([MopSB.BackLat(n) MopSB.OffLat(n)],[MopSB.BackLon(n) MopSB.OffLon(n)],'m-','linewidth',2)
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
Z4=Z;clear Z;
load SolanaPostnourishmentGrid29May24.mat
Z5=Z;clear Z;
load SolanaPostnourishmentGrid07Jun24.mat
Z6=Z;clear Z;
load SolanaPostnourishmentGrid11Jun24.mat
Z7=Z;clear Z;
load SolanaPostnourishmentGrid7Jul24.mat
Z8=Z;clear Z;
%Z4(Z4 >0)=Z4(Z4 > 0)+0.0887;% adding 10cm to AtvMR for datum offset with truck lidar

load SolanaShoreboxMap.mat
idx=find(Z0(:) > -12.5 & Z0(:) < 6 & ~isnan(Z0(:)) & ~isnan(Z1(:)) & ~isnan(Z2(:)) & ~isnan(Z3(:)) & ~isnan(Z4(:)) & ~isnan(Z5(:)) & ~isnan(Z6(:)) & ~isnan(Z7(:)) & ~isnan(Z8(:)));
xg2=repmat(xgimg,1670,1);

% make z datum correction of surveys relative to survey Z0 based on pad elevations
ndx=find(Z0(idx) > 4.4 & xg2(idx) >= 800 & xg2(idx) < 2700);
Z1dz=mean(Z0(idx(ndx))-Z1(idx(ndx)))
Z2dz=mean(Z0(idx(ndx))-Z2(idx(ndx)))
Z3dz=mean(Z0(idx(ndx))-Z3(idx(ndx)))
Z4dz=mean(Z0(idx(ndx))-Z4(idx(ndx)))
Z5dz=mean(Z0(idx(ndx))-Z5(idx(ndx)))
Z6dz=mean(Z0(idx(ndx))-Z6(idx(ndx)))
Z7dz=mean(Z0(idx(ndx))-Z7(idx(ndx)))
Z8dz=mean(Z0(idx(ndx))-Z8(idx(ndx)))

dZ1=Z0*NaN;dZ1(idx)=Z1(idx)-Z0(idx)+Z1dz;
dZ2=Z0*NaN;dZ2(idx)=Z2(idx)-Z0(idx)+Z2dz;
dZ3=Z0*NaN;dZ3(idx)=Z3(idx)-Z0(idx)+Z3dz;
dZ4=Z0*NaN;dZ4(idx)=Z4(idx)-Z0(idx)+Z4dz;
dZ5=Z0*NaN;dZ5(idx)=Z5(idx)-Z0(idx)+Z5dz;
dZ6=Z0*NaN;dZ6(idx)=Z6(idx)-Z0(idx)+Z6dz;
dZ7=Z0*NaN;dZ7(idx)=Z7(idx)-Z0(idx)+Z7dz;
dZ8=Z0*NaN;dZ8(idx)=Z8(idx)-Z0(idx)+Z8dz;
dZL=Z0*NaN;dZL(idx)=Z8(idx)-Z7(idx)-Z7dz+Z8dz;

hold on;
imagesc(xgimg,ygimg,dZL,'AlphaData',~isnan(dZ8));hold on;
set(gca,'clim',[-1 1]);colormap(flipud(polarmap));cb=colorbar;
cb.Label.String='Elevation Change (m)';
 cb.Ticks=[-1:.2:1];
set(gca,'ylim',[-200 850]);
set(gca,'position',[0.1259    0.22500    0.7503    0.8150]);


% [C,h] =contour(xgimg,ygimg,Z0P-0.775,[0 -4 -10],'m:','linewidth',2);
% text(-275,500,{'PreNourishment',' 0, -4, -10m MSL'},'color','m','fontsize',18,'fontweight','bold')
set(gca,'ylim',[-100 250]);
set(gca,'dataaspectratio',[1 .7 1]);hold on;
plot([400 400],[-100 100],'y-','linewidth',5);plot([800 800],[-100 100],'y-','linewidth',5);
plot([2700 2700],[-100 0],'y-','linewidth',5);plot([3000 3000],[-100 0],'y-','linewidth',5);
text(600,-30,{'South of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(2850,-30,{'North of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(1750,-10,{'Pad Reach'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
title('Solana Beach Nourishment Pad Evolution :  21 Jun (ATV Miniranger) vs 11 Jun 2024','fontsize',20)


%%
ax11=axes;
hold on;
image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');hold on;
set(gca,'dataaspectratio',[1 1 1],'xdir','normal','fontsize',16);
set(gca,'ylim',[min(ygimg) max(ygimg)],'xlim',[min(xgimg) max(xgimg)],'xtick',[])
pos=get(gca,'position');
%set(gca,'position',[.5 pos(2) pos(3) pos(4)],'fontsize',14,'linewidth',2)
%xlabel('Alongshore Distance (m)','fontsize',18);
%ylabel('Cross-shore Distance (m)','fontsize',18);

imagesc(xgimg,ygimg,dZ8,'AlphaData',~isnan(dZ5));hold on;
set(gca,'clim',[-1 1]);colormap(flipud(polarmap));cb=colorbar;
cb.Label.String='Elevation Change (m)';
 cb.Ticks=[-1:.2:1];
set(gca,'ylim',[-200 850]);
set(gca,'position',[0.1259    0.4600    0.7503    0.8150]);

% [C,h] =contour(xgimg,ygimg,Z0P-0.775,[0 -4 -10],'m:','linewidth',2);
% text(-275,500,{'PreNourishment',' 0, -4, -10m MSL'},'color','m','fontsize',18,'fontweight','bold')
set(gca,'ylim',[-100 250]);
set(gca,'dataaspectratio',[1 .7 1]);hold on;
plot([400 400],[-100 100],'y-','linewidth',5);plot([800 800],[-100 100],'y-','linewidth',5);
plot([2700 2700],[-100 0],'y-','linewidth',5);plot([3000 3000],[-100 0],'y-','linewidth',5);
text(600,-30,{'South of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(2850,-30,{'North of','Pad'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
text(1750,-10,{'Pad Reach'},'horizontalalign','center','fontsize',20,'color','y',...
    'fontweight','bold','backgroundcolor',[.7 .7 .7]);
title('Solana Beach Nourishment Pad Evolution :  11 Jun (ATV Miniranger) vs 16 Apr 2024','fontsize',20)

%clabel(C,h)
% title({'Solana Beach Nourishment',[datestr(sdate) ' Truck LiDAR Survey'],...
%     'Observed Vol Change above MSL^{**} = 440K m^3 (525K planned) ** new subaqueous toe not included'})
%%
ax2=axes('position',[0.1259    0.3100    0.7503    0.1650]);

p1=plot(xgimg,sum(dZ1,'omitnan'),'-','linewidth',2);hold on;
p2=plot(xgimg,sum(dZ2,'omitnan'),'-','linewidth',2);hold on;
p3=plot(xgimg,sum(dZ3,'omitnan'),'-','linewidth',2);hold on;
p4=plot(xgimg,sum(dZ4,'omitnan'),'-','linewidth',2);hold on;
p5=plot(xgimg,sum(dZ5,'omitnan'),'-','linewidth',2);hold on;
p6=plot(xgimg,sum(dZ6,'omitnan'),'-','linewidth',2);hold on;
p7=plot(xgimg,sum(dZ7,'omitnan'),'-','linewidth',2);hold on;
p8=plot(xgimg,sum(dZ8,'omitnan'),'-','linewidth',2);hold on;
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
SP4=sprintf('%+8i',round(sum(sum(dZ4(:,idx),'omitnan'))));
dSP4=sprintf('(%+8i)',round(sum(sum(dZ4(:,idx),'omitnan')))-round(sum(sum(dZ3(:,idx),'omitnan'))));
SP5=sprintf('%+8i',round(sum(sum(dZ5(:,idx),'omitnan'))));
dSP5=sprintf('(%+8i)',round(sum(sum(dZ5(:,idx),'omitnan')))-round(sum(sum(dZ4(:,idx),'omitnan'))));
SP6=sprintf('%+8i',round(sum(sum(dZ6(:,idx),'omitnan'))));
dSP6=sprintf('(%+8i)',round(sum(sum(dZ6(:,idx),'omitnan')))-round(sum(sum(dZ5(:,idx),'omitnan'))));
SP7=sprintf('%+8i',round(sum(sum(dZ7(:,idx),'omitnan'))));
dSP7=sprintf('(%+8i)',round(sum(sum(dZ7(:,idx),'omitnan')))-round(sum(sum(dZ6(:,idx),'omitnan'))));
SP8=sprintf('%+8i',round(sum(sum(dZ8(:,idx),'omitnan'))));
dSP8=sprintf('(%+8i)',round(sum(sum(dZ8(:,idx),'omitnan')))-round(sum(sum(dZ7(:,idx),'omitnan'))));
% exclude > pad elevations from pad reach change calc
zcut=4.4;% cutoff for AtvMR;% zcut=4.4; % truc lidar
dZ1p=dZ1;dZ1p(Z1 > zcut)=NaN;
dZ2p=dZ2;dZ2p(Z2 > zcut)=NaN;
dZ3p=dZ3;dZ3p(Z3 > zcut)=NaN;
dZ4p=dZ4;dZ4p(Z4 > zcut)=NaN;
dZ5p=dZ5;dZ5p(Z5 > zcut)=NaN;
dZ6p=dZ6;dZ6p(Z6 > zcut)=NaN;
dZ7p=dZ7;dZ7p(Z7 > zcut)=NaN;
dZ8p=dZ7;dZ8p(Z8 > zcut)=NaN;
idx=find(xgimg >= 800 & xgimg < 2700);
P1=sprintf('%+8i',round(sum(sum(dZ1p(:,idx),'omitnan'))));
P2=sprintf('%+8i',round(sum(sum(dZ2p(:,idx),'omitnan'))));
dP2=sprintf('(%+8i)',round(sum(sum(dZ2p(:,idx),'omitnan')))-round(sum(sum(dZ1p(:,idx),'omitnan'))));
P3=sprintf('%+8i',round(sum(sum(dZ3p(:,idx),'omitnan'))));
dP3=sprintf('(%+8i)',round(sum(sum(dZ3p(:,idx),'omitnan')))-round(sum(sum(dZ2p(:,idx),'omitnan'))));
P4=sprintf('%+8i',round(sum(sum(dZ4p(:,idx),'omitnan'))));
dP4=sprintf('(%+8i)',round(sum(sum(dZ4p(:,idx),'omitnan')))-round(sum(sum(dZ3p(:,idx),'omitnan'))));
P5=sprintf('%+8i',round(sum(sum(dZ5p(:,idx),'omitnan'))));
dP5=sprintf('(%+8i)',round(sum(sum(dZ5p(:,idx),'omitnan')))-round(sum(sum(dZ4p(:,idx),'omitnan'))));
P6=sprintf('%+8i',round(sum(sum(dZ6p(:,idx),'omitnan'))));
dP6=sprintf('(%+8i)',round(sum(sum(dZ6p(:,idx),'omitnan')))-round(sum(sum(dZ5p(:,idx),'omitnan'))));
P7=sprintf('%+8i',round(sum(sum(dZ7p(:,idx),'omitnan'))));
dP7=sprintf('(%+8i)',round(sum(sum(dZ7p(:,idx),'omitnan')))-round(sum(sum(dZ6p(:,idx),'omitnan'))));
P8=sprintf('%+8i',round(sum(sum(dZ8p(:,idx),'omitnan'))));
dP8=sprintf('(%+8i)',round(sum(sum(dZ8p(:,idx),'omitnan')))-round(sum(sum(dZ7p(:,idx),'omitnan'))));
idx=find(xgimg >= 2700 & xgimg < 3001);
NP1=sprintf('%+8i',round(sum(sum(dZ1(:,idx),'omitnan'))));
NP2=sprintf('%+8i',round(sum(sum(dZ2(:,idx),'omitnan'))));
dNP2=sprintf('(%+8i)',round(sum(sum(dZ2(:,idx),'omitnan')))-round(sum(sum(dZ1(:,idx),'omitnan'))));
NP3=sprintf('%+8i',round(sum(sum(dZ3(:,idx),'omitnan'))));
dNP3=sprintf('(%+8i)',round(sum(sum(dZ3(:,idx),'omitnan')))-round(sum(sum(dZ2(:,idx),'omitnan'))));
NP4=sprintf('%+8i',round(sum(sum(dZ4(:,idx),'omitnan'))));
dNP4=sprintf('(%+8i)',round(sum(sum(dZ4(:,idx),'omitnan')))-round(sum(sum(dZ3(:,idx),'omitnan'))));
NP5=sprintf('%+8i',round(sum(sum(dZ5(:,idx),'omitnan'))));
dNP5=sprintf('(%+8i)',round(sum(sum(dZ5(:,idx),'omitnan')))-round(sum(sum(dZ4(:,idx),'omitnan'))));
NP6=sprintf('%+8i',round(sum(sum(dZ6(:,idx),'omitnan'))));
dNP6=sprintf('(%+8i)',round(sum(sum(dZ6(:,idx),'omitnan')))-round(sum(sum(dZ5(:,idx),'omitnan'))));
NP7=sprintf('%+8i',round(sum(sum(dZ7(:,idx),'omitnan'))));
dNP7=sprintf('(%+8i)',round(sum(sum(dZ7(:,idx),'omitnan')))-round(sum(sum(dZ6(:,idx),'omitnan'))));
NP8=sprintf('%+8i',round(sum(sum(dZ8(:,idx),'omitnan'))));
dNP8=sprintf('(%+8i)',round(sum(sum(dZ8(:,idx),'omitnan')))-round(sum(sum(dZ7(:,idx),'omitnan'))));


% str = regexprep(str,'(?<!\.\d*)\d{1,3}(?=(\d{3})+\>)','$&,');
% 
% str2=num2str(round(sum(sum(dZ3(:,700:end),'omitnan'))),'%i');
% str2 = regexprep(str2,'(?<!\.\d*)\d{1,3}(?=(\d{3})+\>)','$&,');

% text(500,60,['Total Change ' str ' m^3  (' str2 ' m^3 excluding Dog Beach)'],...
%     'backgroundcolor','w','fontsize',22)
%legend([p2 p1],'16 Apr to 01 May','16 Apr to 10 May','location','northeast')

%% --------- 

ax3=axes('position',[0.1259    0.065    0.7503    0.1650]);
% vertical resolution for summed changes
Res=0.1;
% round survey x,y to desired resolution
%zr=Res*round(z1da/Res); % round to Res meter spatial resolution
zr=Res*round(Z0(:)/Res);

% keep changes for the pad reach only
dz=dZ1(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b1=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 26 Apr[' num2str(-Z1dz*100,'%+4.1f') ']' SP1 '          ' P1 '             ' NP1]);hold on;

dz=dZ2(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b2=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 01 May[' num2str(-Z2dz*100,'%+4.1f') ']' SP2 dSP2 P2 dP2 '  ' NP2 dNP2]);hold on;

dz=dZ3(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b3=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 10 May[' num2str(-Z3dz*100,'%+4.1f') ']' SP3 dSP3 ' ' P3 dP3 ' ' NP3 dNP3]);hold on;

dz=dZ4(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b4=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 24 May[' num2str(-Z4dz*100,'%+4.1f') ']' SP4 dSP4 ' ' P4 dP4 ' ' NP4 dNP4]);hold on;

dz=dZ5(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b5=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 29 May[' num2str(-Z5dz*100,'%+4.1f') ']' SP5 dSP5 ' ' P5 dP5 ' ' NP5 dNP5]);hold on;


dz=dZ6(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b6=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 07 Jun[' num2str(-Z6dz*100,'%+4.1f') ']' SP6 dSP6 ' ' P6 dP6 ' ' NP6 dNP6]);hold on;


dz=dZ7(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b7=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 11 Jun[' num2str(-Z7dz*100,'%+4.1f') ']' SP7 dSP7 ' ' P7 dP7 ' ' NP7 dNP7]);hold on;


dz=dZ8(:);dz(X(:) < 800)=NaN;dz(X(:) > 2700)=NaN;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
b8=plot(uz-0.774,dztot/1000,'-','linewidth',2,'displayname',...
    ['16 Apr to 27 Jun[' num2str(-Z8dz*100,'%+4.1f') ']' SP8 dSP8 ' ' P8 dP8 ' ' NP8 dNP8]);hold on;

ytickformat('%g K')
set(gca,'xlim',[-12 5],'ylim',[-.4 .4])
set(gca,'fontsize',16);grid on;%set(gca,'ytick',0:100:400,'xlim',[xgimg(1) xgimg(end)])
xlabel('16 Apr Elevation (m, MSL)','fontsize',18);
ylabel({'Vol Change (m^{3})'})
set(gca,'ytick',-0.4:.2:.4)
text(-11.5,-.25,{'Pad Reach Volume Change','   vs. 16 Apr Elevations'},...
    'backgroundcolor','w','fontsize',16)
hold on;
text(16.2,.52,'Reach Volume Change (m^{3})','fontsize',22)
text(6.,-.52,'( )=Change from Previous Survey ;** Pad Terrace Excluded from Vol Change Estimate','fontsize',16)
%text(6.,-.55,'** Pad Terrace Excluded from Vol Change Estimate','fontsize',18)
text(6.,-.64,'* Z Datum Offset Based on Mean Pad Terrace Elevations','fontsize',16)


plot([-.83 -.83],[-.4 .4],'k--','linewidth',2)
plot([.79 .79],[-.4 .4],'k--','linewidth',2)
plot([3.5 3.5],[-.4 .4],'k:','linewidth',2)
plot([4 4],[-.4 .4],'k:','linewidth',2)
text([-.88],[.4],'MLLW','fontsize',14,'HorizontalAlignment','right')
text([.86],[.4],'MHHW','fontsize',14,'HorizontalAlignment','left')
text([3],[.54],{'   Pad','Terrace'},'color','m','fontsize',14,'HorizontalAlignment','left')

set(gcf,'color','w')
set(gcf,'inverthardcopy','off')
set(gcf,'PaperPositionMode','auto');
lg=legend([b1 b2 b3 b4 b5 b6 b7 b8],'location','eastoutside','fontsize',14,'fontname','Monospaced');
title(lg,'Date [z offset(cm)*]         South         Pad Reach^{**}          North ')
print(gcf,'-dpng','-r300','-loose','SolanaNourishment27JunVs16Apr24.png')
