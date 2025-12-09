load MopLidarVolumesMSLib.mat

% make a mop number to go with every survey date
for n=1:size(LidarVols,2)
    LidarVols(n).Mop=repmat(LidarVols(n).Mop,size(LidarVols(n).Vdatenums));
end
%%
% 
figure;plot3(vertcat([LidarVols.Vdatenums]),vertcat([LidarVols.Mop]),vertcat([LidarVols.VolMin]),'.')
% MobileLidarVols=vertcat(LidarVols.VolMin);

%% plot min surface
figure('position',[246   262   871   387]);

y=vertcat([LidarVols.Mop]);
z=vertcat([LidarVols.VolMin]/100);

zqc=z;

x=[LidarVols.Vdatenums];
for ym=min(y):max(y)
    idx=find(y == ym);
    mdx=find(y == ym & x > datenum(2019,1,1) & x < datenum(2024,1,1));
    mdx=idx; % ib
    bdx=find(isoutlier(z(idx)));
    zqc(idx(bdx))=NaN;
end

%z=zqc;
z=vertcat([LidarVols.VolMin]/100);
x=[LidarVols.Vdatetimes];

% plot colored points
minz=0;%min(z);%quantile(z,.05);
maxz=200;%max(z);% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(x(idx), y(idx), 15, cm(ceil(zscaled(idx)),:), 'filled','s');
cb=colorbar;cb.Label.String=' MSL Subaerial Mobile Volume (m^{3}/m-shoreline) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
set(gca,'color','w','xlim',[datetime(2001,1,1) datetime(2025,1,1)]);grid on;
set(gca,'xtick',datetime(2001:2025,1,1))
ylabel('MOP')
makepng('IBLidarSubaerialVolumes.png')
% set(gca,'xlim',[BeachAreaLon-0.0012 BeachAreaLon+0.0012]);
% set(gca,'ylim',[BeachAreaLat-0.0012 BeachAreaLat+0.0012]);
% set(gca,'xtick',[],'ytick',[])

%% 
figure
% z=vertcat([LidarVols.VolMin]/100);
% zqc=z;
zanom=z*NaN;

x=[LidarVols.Vdatenums];
n=0;
for ym=min(y):max(y)
    idx=find(y == ym);
    mdx=find(y == ym & x > datenum(2019,1,1) & x < datenum(2024,1,1)); % north county
    mdx=idx;  % imperial beach
    bdx=find(isoutlier(z(idx)));
    zqc(idx(bdx))=NaN;
    n=n+1;Vmean(n)=mean(z(mdx),'omitnan');Vstd(n)=std(z(mdx),'omitnan');
    zanom(idx)=z(idx)-mean(z(mdx),'omitnan');
end

x=[LidarVols.Vdatetimes];
y=vertcat([LidarVols.Mop]);
z=zanom;
z(z < 0)=NaN;
% plot colored points
minz=0;%quantile(z,.05);
maxz=75;% zmax=quantile(z,.95);
zrange=maxz-minz; % set max slope for coloring
cm =jet(64);%flipud(polarmap(64));
zscaled = 1+64*(z-minz)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
scp=scatter(x(idx), y(idx), 8, cm(ceil(zscaled(idx)),:), 'filled');
cb=colorbar;cb.Label.String=' Subaerial Mobile Volume (m^{3}/m-shoreline) ';
set(gca,'clim',[minz maxz],'fontsize',16);colormap(cm);
set(gca,'color','w','xlim',[datetime(2000,1,1) datetime(2025,1,1)])


figure
%errorbar(min(y):max(y),Vmean,Vstd)
Ax=axes;
Ax.Layer = 'top';
Ax.GridAlpha = 0.5;
p(2)=fill([min(y):max(y) flip(min(y):max(y))], [Vmean-Vstd flip(Vmean+Vstd)],[.88 .88 .88],'facealpha',1,'EdgeColor','none',....
    'displayname','Std Dev All Transects');
hold on;
plot(min(y):max(y),Vmean,'k-','linewidth',2)
ylabel('MSL Subaerial Mobile Volume Mean & Std Dev (m^{3}/m-shoreline)');
xlabel('MOP');
grid on
makepng('IBmeanLidarSubaerialVolumes.png')