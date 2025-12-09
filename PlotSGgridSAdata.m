% shoals 9/4/14 survey
%iyr=2001;imn=11;idy=28;
load M00582SG.mat
ndx=find(strcmpi({SG.Source},'Gps') == 1 );
imop=310;
SAT=[];
for MopNumber=580:584%570:595
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat' ],'SA');    
idx=find(strcmp({SA.File},SG(imop).File) == 1);
SAT=[SAT SA(idx)];
end

figure;
%SA=SAT;
clf;ScatterPlotBeachBarUTM(vertcat(SAT.X),vertcat(SAT.Y),vertcat(SAT.Z),'3d')
grid on;
zlabel('Elevation (m, NAVD88)')
BeachBarColorbar
set(gca,'fontsize',14)
%title("Connor's Most Excellent First Beach Survey")
title(datestr(SG(imop).Datenum));
view(2)
for MopNumber=580:584
    PlotLabelMopTransectUTM(MopNumber,'3d','k','ShadeOff')
end

% %----------
% 
% Xr=vertcat(SAT.X);
% Yr=vertcat(SAT.Y);
% Zr=vertcat(SAT.Z);
% 
%      fprintf('Gridding %i points...\n',length(Xr))
%      XGutm=[];YGutm=[];Vzg=[];
%     
% %      min(Yr)
% %      max(Yr)
% % bound the data gaps with NaNs to avoid gridding them
%      MaxGap=250;
%      %Xr(Yr > 3650000)=[];Zr(Yr > 3650000)=[];Yr(Yr > 3650000)=[];
%      [x,y,z]=addNoDataAreaPoints(Xr,Yr,Zr,MaxGap);     
%     
% % Grid the elevation survey using Delaunay tesselation 
%     zg=griddata(double(x),double(y),double(z),...
%         double(min(x):max(x)),double(min(y):max(y))');
%     
%     [YG,XG]=find(~isnan(zg)); % find valid grid data x,y vectors
%      Vzg=zg(~isnan(zg(:))); % valid zg data vector
%      XGutm=min(x)-1+XG; % adjust back to utm coords
%      YGutm=min(y)-1+YG;
%      
% figure;
% ScatterPlotBeachBarUTM(XGutm,YGutm,Vzg,'3d');
% grid on;
% zlabel('Elevation (m, NAVD88)')
% BeachColorbar
% set(gca,'fontsize',14)