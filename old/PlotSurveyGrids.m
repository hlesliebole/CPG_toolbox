%PlotSurveyGrids.m
%close all

%  SUPER slow with big struct arrays
figure
for n=1:size(SG,2)
    p(n)=scatter3(SG(n).X,SG(n).Y,SG(n).Z,75,'.');
    p(n).MarkerFaceColor = p(n).CData;
    p(n).MarkerFaceAlpha = 1;
    p(n).MarkerEdgeAlpha = 1;
    hold on;
end
set(gca,'dataaspectratio',[1 1 .05]);
xl=get(gca,'xlim');yl=get(gca,'ylim');
zd=0.744; % msl
fill3([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],...
    [zd zd zd zd zd],[.1 .1 .1],...
        'edgecolor','none','facealpha',.2)
    zd=1.566; % mhhw
fill3([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],...
    [zd zd zd zd zd],[.1 .1 .1],...
        'edgecolor','none','facealpha',.2)
legend(p,datestr(vertcat(SG.Datenum)),'location','northwest')
    