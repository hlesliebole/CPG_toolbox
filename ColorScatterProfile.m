function [ScatterPlot,ColorBarPlot]=ColorScatterProfile(x,y,z)

[z,i]=sort(z,'descend');
x=x(i);
y=y(i);

zmin=quantile(z,.05);zmax=quantile(z,.95);

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

cb=colorbar;cb.Label.String='Dist From Mop Line (m)';
set(gca,'clim',[zmin zmax]);
set(gca,'xlim',[min(x) max(x)]);
%set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end