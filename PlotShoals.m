% shoals 9/4/14 survey
%iyr=2014;imn=9;idy=4;
iyr=2009;imn=10;idy=12;
figure('position', [279    57   800   734]);

SAT=[];
for MopNumber=585:595%560:580%570:595
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat' ],'SA');    
idx=find([SA.Datenum] == datenum(iyr,imn,idy));
SAT=[SAT SA(idx)];
end

%SA=SAT;
clf;ScatterPlotBeachUTM(vertcat(SAT.X),vertcat(SAT.Y),vertcat(SAT.Z),'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
%title("Shoals 9/4/2014 Survey : Inlet Mops 585-595")
title("Shoals 10/12/2009 Survey : Inlet Mops 585-595")

view(2)
% for MopNumber=560:580
%     PlotLabelMopTransectUTM(MopNumber,'3d','k','ShadeOff')
% end