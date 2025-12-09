
load M00550SA.mat
SAT=SA(end);
load M00551SA.mat
SAT=[SAT SA(end)];
load M00552SA.mat
SAT=[SAT SA(end)];
load M00553SA.mat
SAT=[SAT SA(end)];
load M00554SA.mat
SAT=[SAT SA(end)];
load M00555SA.mat
SAT=[SAT SA(end)];

SA=SAT;
figure('position',[191   261   846   491]);
ScatterPlotBeachUTM(vertcat(SA.X),vertcat(SA.Y),vertcat(SA.Z),'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
title([ datestr(SA(end).Datenum,'mm/dd/yyyy') " Blacks North Wheelie Survey"])
makepng('LatestWheelieBN.png')