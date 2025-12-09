load M00690SM.mat
figure;
hold on;plot(SM(5).X1D,SM(5).X1D*0+0.774,'b-','linewidth',2);
text(-50,1.2,'MSL')
set(gca,'xdir','reverse','fontsize',12)
grid on;
xlabel('Cross-shore Distance (m)');ylabel('Depth (m, NAVD88)');
title([{'MOP 690 , San Elijo State Beach'},...
    {'USACE SHOALS Airborne Bathymetric Survey : ' datestr(SM(5).Datenum)}])

msl=0.774;
x1=SM(5).X1D(1);x2=SM(5).X1D(end);
y1=msl-1;y2=msl-3; % seagrass zone
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'g','facealpha',.2);

y1=msl-3;y2=msl-8; % algae zone
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'y','facealpha',.2);

y1=msl-8;y2=msl-18; % kelp zone
fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r','facealpha',.2);

set(gca,'xlim',[x1 x2],'ylim',[-18 5]);

text(50,msl-2,'SEAGRASS ZONE','verticalalign','middle','fontweight','bold');
text(50,msl-5.5,'ALGAE ZONE','verticalalign','middle','fontweight','bold');
text(50,msl-13,'KELP ZONE','verticalalign','middle','fontweight','bold');

plot(SM(5).X1D,SM(5).Z1Dmean,'k-','linewidth',2)
print(gcf,'-dpng','-r300','-loose','SanElijoSBplantZones.png');

