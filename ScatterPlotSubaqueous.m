function ScatterPlotSubaqueousUTM(x,y,z,zmin,PlotMethod)

% makes jet colormap scatter plot for depths between zmin and 0.774 (MSL)

x(z > 0.774)=[];y(z > 0.774)=[];z(z > 0.774)=[];
%load BeachColorMap
BeachColorMap=jet(64);
BeachColorBarLims=[zmin 0.774];

zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    (BeachColorBarLims(2)-BeachColorBarLims(1));
zscaled(zscaled < 1)=1;
zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);
zscaled(isnan(z))=1;

%cn = ceil(max(zscaled));                                       
cm = BeachColorMap;

if contains(PlotMethod,'3d','IgnoreCase',true) % 3D plot
%scp=scatter3(x, y, z, 10, cm(ceil(zscaled),:), 'filled');
scp=scatter3(x, y, z, 1, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;

else
%scp=scatter(x, y, 10, cm(ceil(zscaled),:), 'filled');
scp=scatter(x, y, 1, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;
end

%axis equal
set(gca,'fontsize',14);grid;box on;
%axis equal;
% if(min(z) > -4)
%  set(gca,'dataaspectratio',[1 1 .03]);
% else
%  set(gca,'dataaspectratio',[1 1 .05]);
% end

y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('Latitude');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('Longitude');

colormap(jet(64));
set(gca,'clim',[zmin 0.774]);
cb=colorbar;
cb.Label.String='Elevation (m, NAVD88)';
cb.FontSize=14;

end
