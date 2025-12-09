% example code to make some plots of MHW beach widths
%  from the MopShoreline table array

load MopShorelineV2.mat

x=[];y=[];zw=[];za=[];

mnum=MopShoreline.MopNum;
for m=3:920 % border to oceanside
    idx=find(mnum == m); % find this mop's rows in table
    d=MopShoreline.DateNum(idx); %survey datnums
%     xMSL=MopShoreline.MslDist(idx); 
%     TF=isoutlier(xMSL);xMSL(TF)=NaN; 
    xMHW=MopShoreline.MhwDist(idx); % mhw widths
    TF=isoutlier(xMHW);xMHW(TF)=NaN; % remove outliers
    % get season-weighted mean width of all surveys
    [xMhwMean,qmean,mmean,ymmean]=SeasonWeightedMean(d,xMHW);     
    aMHW=xMHW-xMhwMean; % width anomaly
    % add to vectors of data for all the mops
    x=[x d'];y=[y mnum(idx)'];za=[za aMHW'];zw=[zw xMHW'];
end

%% Fig 1 Plot most recent width at each Mop relative to a gray cloud
%        of past widths
figure('position',[ 7  101 1287   700]);
p1=plot(y,zw,'.','color',[.8 .8 .8]);hold on
% find any widths after Jan 6, 2023
idx=find(x > datenum(2023,1,6));
p2=plot(y(idx),zw(idx),'r.','markersize',10);
set(gca,'xlim',[0 925],'ylim',[-50 250],'fontsize',14);
grid on
ylabel('MHW Beach Width (m)');xlabel('Mop Number');
legend([p1 p2],'Past Observed Widths','Widths Since 6 Jan 2023')

%% Fig 2 Plot most recent width anomaly at each Mop relative to a gray cloud
%        of past width anomalies
figure('position',[ 27  95 1287   700]);
p1=plot(y,za,'.','color',[.8 .8 .8]);hold on
% find any widths after Jan 6, 2023
idx=find(x > datenum(2023,1,6));
p2=plot(y(idx),za(idx),'r.','markersize',10);
plot([0 925],[0 0],'k-');
set(gca,'xlim',[0 925],'ylim',[-60 60],'fontsize',14);
grid on
ylabel('MHW Beach Width Anomaly (m)');xlabel('Mop Number');
legend([p1 p2],'Past Observed Width Anomaly','Width Anomalies Since 6 Jan 2023')

%% for Figs 3,4 sort widths and anomalies in descending order to make sure red narrow
%  and erosion anomalies are seen
[za,is]=sort(za,'descend');xa=x(is);ya=y(is);
[zw,is]=sort(zw,'descend');x=x(is);y=y(is);
           
%% Fig 3 make 2d mhw width plot
figure('position',[ 67          81        1287         700]);

cm =flipud(jet(64)); % jet colormap
[scp,cb]=ColorScatter(x,y,zw,cm);

cb.Label.String='MHW Beach Width(m)';
cb.FontSize=14;
set(gca,'xlim',[datenum(1997,1,0) datenum(2024,1,0)])
set(gca,'ylim',[0 925],'fontsize',14);
datetick;view(2);grid on
xlabel('Survey Date');ylabel('Mop Number');

%% Fig 4 make 2d anomaly plot
figure('position',[ 111          36        1287         700]);

% make polarmap scatter plot of data points
cm =flipud(polarmap(64)); % polar colormap
[scp,cb]=ColorScatterPolar(xa,ya,za,cm);

cb.Label.String='MHW Beach Width Anomaly (m)';
cb.FontSize=14;
set(gca,'xlim',[datenum(1997,1,0) datenum(2024,1,0)])
set(gca,'ylim',[0 925],'fontsize',14);
datetick;view(2);grid on
xlabel('Survey Date');ylabel('Mop Number');
%set(gca,'view',[51.6348   74.0140]); % alternative 3d view

%% ------------------------------------------------------------------

%%
function [swmean,qmeans,mmeans,ymmeans]=SeasonWeightedMean(d,xMHW)

dtime=datetime(d,'convertfrom','datenum');
y=year(dtime);
mn=month(dtime);
dy=day(dtime);

% first reduce to individual year-month means
uy=unique(y); %unique years
ymmeans(1:numel(uy(1):uy(end)),1:12)=NaN;
for ny=uy'
    for m=unique(mn(y == ny))'
        ymmeans(ny-uy(1)+1,m)=mean(xMHW(y == ny & mn == m),'omitnan');
    end
end
        
% then reduce to month means
mmeans(1:12)=NaN;
for m=unique(mn)'
    mmeans(m)=mean(ymmeans(:,m),'omitnan');
end

% now quarterly means
qmeans(1:4)=NaN;
for q=1:4
    qmeans(q)=mean(mmeans((q-1)*3+1:(q-1)*3+3),'omitnan');
end

% finally season weighted global mean
swmean=mean(qmeans,'omitnan');

end

%%
function [ScatterPlot,ColorBarPlot]=ColorScatterPolar(x,y,z,cm)

%cm =flipud(polarmap(64));

zmax=max([abs(quantile(z,.05)) quantile(z,.95)]);
zmin=-zmax;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
scp=scatter3(x(idx), y(idx), z(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-10,80)
colormap(flipud(polarmap(64)))


set(gca,'clim',[zmin zmax]);
set(gca,'xlim',[min(x) max(x)]);
set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);
cb=colorbar;

ScatterPlot=scp;
ColorBarPlot=cb;

end
%%
function [ScatterPlot,ColorBarPlot]=ColorScatter(x,y,z,cm)

%cm =flipud(polarmap(64));

zmin=quantile(z,.02);
zmax=quantile(z,.98);
%zmax=max([abs(quantile(z,.05)) quantile(z,.95)]);
%zmin=-zmax;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
scp=scatter3(x(idx), y(idx), z(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-10,80)
colormap(cm)


set(gca,'clim',[zmin zmax]);
set(gca,'xlim',[min(x) max(x)]);
set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);
cb=colorbar;

ScatterPlot=scp;
ColorBarPlot=cb;

end