clearvars
load MopSubaqueousVolumesMSL.mat

% make a mop number to go with every survey date
for n=1:size(TransectVols,2)
    TransectVols(n).Mop=repmat(TransectVols(n).Mop,size(TransectVols(n).Vdatenums));
end
%%
% 
%figure;plot3(vertcat([TransectVols.Vdatenums]),vertcat([TransectVols.Mop]),vertcat([TransectVols.VolTot]),'.')
% MobileTransectVols=vertcat(TransectVols.VolTot);

%% plot min surface
figure('position',[246   262   871   387]);

y=vertcat([TransectVols.Mop]);
z=vertcat([TransectVols.VolTot]/100);

zqc=z;

x=[TransectVols.Vdatenums];
for ym=min(y):max(y)
    idx=find(y == ym);
    mdx=find(y == ym & x > datenum(2003,1,1) & x < datenum(2024,1,1));
    %mdx=idx; % ib
    bdx=find(isoutlier(z(idx)));
    zqc(idx(bdx))=NaN;
end

%z=zqc;
z=vertcat([TransectVols.VolTot]/100);
x=[TransectVols.Vdatetimes];

% plot colored points
minz=100;%min(z);%quantile(z,.05);
maxz=600;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(x(idx), y(idx), 15, cm(ceil(zscaled(idx)),:), 'filled','s');
cb=colorbar;cb.Label.String=' MSL Subaqueous Mobile Volume (m^{3}/m-shoreline) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
set(gca,'color','w','xlim',[datetime(1997,1,1) datetime(2025,1,1)]);grid on;
set(gca,'xtick',datetime(1997:2025,1,1))
ylabel('MOP')

%makepng('LidarSubaerialVolumes.png')
% set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
% set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
% set(gca,'xtick',[],'ytick',[])

figure('position',[ 37         292        1392         429])
scp=scatter(y(idx), x(idx), 15, cm(ceil(zscaled(idx)),:), 'filled','s');
cb=colorbar;cb.Label.String=' MSL Subaqueous Mobile Volume (m^{3}/m-shoreline) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
set(gca,'color','w','ylim',[datetime(1997,1,1) datetime(2025,1,1)]);grid on;
set(gca,'ydir','reverse','ytick',datetime(1997:2025,1,1))
xlabel('MOP')
set(gca,'xlim',[518 770])
set(gca,'color','k')
set(gcf,'inverthardcopy','off')
xlabel('       South Torrey                   North Torrey       S Del Mar   N Del Mar     Solana Bch        Cardiff             Encinitas                      Leucadia           ')
%makepng('SioCanyon2BatiquitosSubaerialVolume.png')

% %% 
% figure
% % z=vertcat([TransectVols.VolTot]/100);
% %zqc=z;
% zanom=z*NaN;
% 
% x=[TransectVols.Vdatenums];
% n=0;
% for ym=min(y):max(y)
%     idx=find(y == ym);
%     mdx=find(y == ym & x > datenum(2019,1,1) & x < datenum(2024,1,1)); % north county
%     %mdx=idx;  % imperial beach
%     bdx=find(isoutlier(z(idx)));
%     zqc(idx(bdx))=NaN;
%     n=n+1;Vmean(n)=mean(zqc(mdx),'omitnan');Vstd(n)=std(zqc(mdx),'omitnan');
%     zanom(idx)=zqc(idx)-mean(z(mdx),'omitnan');
% end
% 
% x=[TransectVols.Vdatetimes];
% y=vertcat([TransectVols.Mop]);
% z=zanom;
% z(z < 0)=NaN;
% % plot colored points
% minz=0;%quantile(z,.05);
% maxz=75;% zmax=quantile(z,.95);
% zrange=maxz-minz; % set max slope for coloring
% cm =jet(64);%flipud(polarmap(64));
% zscaled = 1+64*(z-minz)/zrange;
% zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
% idx=find(~isnan(zscaled)); % non NaN points 
% scp=scatter(x(idx), y(idx), 8, cm(ceil(zscaled(idx)),:), 'filled');
% cb=colorbar;cb.Label.String=' Subaerial Mobile Volume (m^{3}/m-shoreline) ';
% set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
% set(gca,'color','w','xlim',[datetime(2000,1,1) datetime(2025,1,1)])

%%
x=[TransectVols.Vdatenums];
z=vertcat([TransectVols.VolTot]/100);
zqc=z;
n=0;
for ym=min(y):max(y)
    idx=find(y == ym);
    bdx=find(isoutlier(z(idx)));
    zqc(idx(bdx))=NaN;
    mdx=find(y == ym & x > datenum(2019,1,1) & x < datenum(2025,1,1)); % north county
    %mdx=idx;  % imperial beach
    n=n+1;
    Vmean(n)=mean(zqc(mdx),'omitnan');
    Vstd(n)=std(zqc(mdx),'omitnan');
    %zanom(idx)=zqc(idx)-mean(z(mdx),'omitnan');
end


figure('position',[ 4         127        1392         644])
Ax=axes;
Ax.Layer = 'top';
Ax.GridAlpha = 0.5;
Vmean=movmean(Vmean,3,'omitnan');
Vstd=movmean(Vstd,3,'omitnan');
Vmean=Vmean(1:1:end);
Vstd=Vstd(1:1:end);

idx=find(~isnan(Vmean));
yy=min(y):1:max(y);yy=yy(idx);
p(2)=fill([yy flip(yy)], [Vmean(idx)-Vstd(idx) flip(Vmean(idx)+Vstd(idx))],[.88 .88 .88],'facealpha',1,'EdgeColor','none',....
    'displayname',['2003-2023 Std Dev ( ' char(177) '0.56M m^{3} )']);
hold on;
%plot(min(y):max(y),Vmean,'k-','linewidth',3)
% plot spring 1998
% xs=sort(unique(x));
% idx=find(x == xs(2));
% ys=min(y):max(y);
% zs=ys*NaN;
% zs(round(y(idx))-min(y)+1)=zqc(idx);
% tv=sum(zs,'omitnan')*100;
% %p(3)=plot(ys,zs,'r.-','linewidth',2,'displayname','NASA ATM April 1998');
hold on;

p(1)=plot(min(y):1:max(y),Vmean,'k-','linewidth',3,'displayname','2003-2023 Mean ( 1.25M m^{3} Total )');

%errorbar(min(y):max(y),Vmean,Vstd)
title('MSL Subaqueous Mobile Volume 2003-2023 Mean & Std Dev (m^{3}/m-shoreline)');
ylabel('MSL Subaqueous Mobile Volume (m^{3}/m-shoreline)');
xlabel('MOP');
set(gca,'xlim',[518 770],'ylim',[0 500],'fontsize',16)
%set(gca,'inverthardcopy','off')
xlabel('       South Torrey                   North Torrey       S Del Mar   N Del Mar     Solana Bch        Cardiff             Encinitas                      Leucadia           ')
grid on
%legend(p)

%makepng('SioCanyon2BatiquitosSubaerialVolumeStats.png')