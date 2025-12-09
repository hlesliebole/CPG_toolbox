function [ScatterPlot,ColorBarPlot]=ColorScatter2dDeepElev(x,y,z)

%zmin=quantile(z,.05);zmax=quantile(z,.95);
zmin=-16;zmax=-9;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
%cm =jet(64);
cm =demcmap([-12 3]);
cm(1:51,:)=jet(51);
cm=flipud(jet(64));

scp=scatter(x(idx), y(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
%view(-10,80)
%colormap(jet(64))
%demcmap([-12 3])
colormap(cm)
cb=colorbar;cb.Label.String='Elevation (m, navd88)';
set(gca,'clim',[zmin zmax]);
% set(gca,'xlim',[min(x) max(x)]);
% set(gca,'ylim',[min(y) max(y)]);
%set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end