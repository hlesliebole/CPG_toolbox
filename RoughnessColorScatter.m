function [ScatterPlot,ColorBarPlot]=RoughnessColorScatter(x,y,d,z)

%zmin=quantile(z,.05);zmax=quantile(z,.95);
zmin=0;zmax=0.5;

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
cm =jet(64);

scp=scatter3(x(idx), y(idx), d(idx), 5, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-25,30)
colormap(jet(64))

cb=colorbar;cb.Label.String='2D Spatial Standard Deviation (Roughness)';
set(gca,'clim',[zmin zmax]);zlabel('Elev (m,NAVD88)');
set(gca,'xlim',[min(x) max(x)]);xlabel('easting (m)');
set(gca,'zlim',[min(d) max(d)]);ylabel('northing (m)');
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end