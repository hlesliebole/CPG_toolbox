
load M00578SA.mat
SAT=SA(end);
load M00579SA.mat
SAT=[SAT SA(end)];
load M00580SA.mat
SAT=[SAT SA(end)];
load M00581SA.mat
SAT=[SAT SA(end)];
load M00582SA.mat
SAT=[SAT SA(end)];
load M00583SA.mat
SAT=[SAT SA(end)];
load M00584SA.mat
SAT=[SAT SA(end)];
SA=SAT;
figure('position',[191   261   846   491]);
ScatterPlotBeachUTM(vertcat(SA.X),vertcat(SA.Y),vertcat(SA.Z),'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
title([ datestr(SA(end).Datenum,'mm/dd/yyyy') " Wheelie Survey"])
makepng('LatestWheelie.png')