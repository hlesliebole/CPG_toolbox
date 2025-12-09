%f1='/Volumes/group/ig8_wheel/20220215_00525_00555_BN_wheel/filtered_20220215.llnezts.navd88';
f1='/Volumes/group/ig8_wheel/20220215_00525_00555_BN2post_wheel/filtered_20220215.llnezts.navd88';

load M00550SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=SA(idx);
load M00551SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=[SAT SA(idx)];
load M00552SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=[SAT SA(idx)];
load M00553SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=[SAT SA(idx)];
load M00554SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=[SAT SA(idx)];
load M00555SA.mat
idx=find(strcmpi({SA.File},f1) ==1);
SAT=[SAT SA(idx)];

SA=SAT;
figure('position',[191   261   846   491]);
ScatterPlotBeachUTM(vertcat(SA.X),vertcat(SA.Y),vertcat(SA.Z),'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
title([ datestr(SA(1).Datenum,'mm/dd/yyyy') '   Rob Repaired'])
makepng('Compare2BN.png')