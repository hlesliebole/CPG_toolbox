function ScatterPlotBeachLatLon(y,x,z,PlotMethod)

load BeachColorMap

zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    (BeachColorBarLims(2)-BeachColorBarLims(1));
zscaled(zscaled < 1)=1;
zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);
zscaled(isnan(z))=1;

%cn = ceil(max(zscaled));                                       
cm = BeachColorMap;

if contains(PlotMethod,'3d','IgnoreCase',true) % 3D plot
scp=scatter3(x, y, z, 10, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;

else
scp=scatter(x, y, 10, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;
end

%BeachColorbar
%axis equal
%set(gca,'fontsize',12);grid;box on;
%axis equal;
% if(min(z) > -4)
%  set(gca,'dataaspectratio',[1 1 .03]);
% else
%  set(gca,'dataaspectratio',[1 1 .05]);
% end

% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');

end
