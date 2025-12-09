dfile='20220216_00581_00589_TorreyCobble_RTKdrone.xyz';

%d=dlmread(dfile,'',3,1);
d=load(dfile);x=d(:,1);y=d(:,2);z=d(:,3);
x(z<-2)=[];y(z<-2)=[];z(z<-2)=[];

figure('position',[191   261   846   491]);
ScatterPlotBeachUTM(x,y,z,'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
title([ "20220216 Drone Survey"])
makepng('Drone.png')