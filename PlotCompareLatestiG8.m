
load M00578SA.mat
SAT=SA(end-1);
SAT2=SA(end);
load M00579SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
load M00580SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
load M00581SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
load M00582SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
load M00583SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
load M00584SA.mat
SAT=[SAT SA(end-1)];
SAT2=[SAT2 SA(end)];
SA=SAT;
SA2=SAT2;
figure('position',[89         203        1099         519]);
%ScatterPlotBeachUTM(vertcat(SA.X),vertcat(SA.Y),vertcat(SA.Z),'3d')
p1=plot3(vertcat(SA.X),vertcat(SA.Y),vertcat(SA.Z),'k.');
hold on;
p2=plot3(vertcat(SA2.X),vertcat(SA2.Y),vertcat(SA2.Z),'r.');
grid on;
zlabel('Elevation (m, NAVD88)')
%BeachColorbar
set(gca,'fontsize',14)
set(gca,'zlim',[-0.5 4])
%title("Connor's Most Excellent First Beach Survey")
title([ datestr(SA(end).Datenum,'mm/dd/yyyy') " Wheelie Survey"])
legend([p1 p2],'Cart','Wheel');
makepng('CompLatestWheelie.png')