% estimate the volume of the cardiff dune

% relevant surveys
SurveyDatenum3=datenum(2018,8,9);SurveySource3='Trk';% pre-dune but some beach nourishing?
SurveyDatenum2=datenum(2016,3,3);SurveySource2='Gps';% post el nino 
%SurveyDatenum3=datenum(2016,3,22);SurveySource3='KMair';% post el nino 
SurveyDatenum3=datenum(1998,4,8);SurveySource3='UTAir';% post el nino
SurveyDatenum1=datenum(2015,10,6);SurveySource1='KMair';% Melville airborne lidar prior to el nino
%SurveyDatenum2=datenum(2016,9,29);SurveySource2='Gps';% post-el nino fall jumbo
SurveyDatenum4=datenum(2019,9,12);SurveySource4='Trk';% post-dune 
SurveyDatenum5=datenum(2010,9,23);SurveySource5='CCC';% post- el nino
%SurveyDatenum5=datenum(2021,10,11);SurveySource5='Trk';% last fall 
SurveyDatenum6=datenum(2023,1,9);SurveySource6='Trk';% most recent

% plot profiles for those 4 dates at 676

figure('position',[440   304   591   493]);
MopNumber=670;
matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat' ];
fprintf('Loading %s\n',matfile);
load(matfile,'SM');

% plot road
SurvNum=find([SM.Datenum] == SurveyDatenum5 & strcmp({SM.Source},SurveySource5));
n=find(SM(SurvNum).X1D < 0);
p5=plot(SM(SurvNum).X1D(n),SM(SurvNum).Z1Dtransect(n),'-','linewidth',2,'color',[.8 .8 .8],...
    'DisplayName',[datestr(SurveyDatenum5) ' ' SurveySource5]);
hold on;
text(-17,5.2,'ROAD','color',[.4 .4 .4],'fontsize',12,'fontweight','bold')

SurvNum=find([SM.Datenum] == SurveyDatenum1 & strcmp({SM.Source},SurveySource1));
p1=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',3,...
    'DisplayName',[datestr(SurveyDatenum1) ' ' SurveySource1]);
hold on;
SurvNum=find([SM.Datenum] == SurveyDatenum2 & strcmp({SM.Source},SurveySource2));
p2=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
    'DisplayName',[datestr(SurveyDatenum2) ' ' SurveySource2]);
% 
SurvNum=find([SM.Datenum] == SurveyDatenum3 & strcmp({SM.Source},SurveySource3));
p3=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
    'DisplayName',[datestr(SurveyDatenum3) ' ' SurveySource3]);
% 
% 
SurvNum=find([SM.Datenum] == SurveyDatenum4 & strcmp({SM.Source},SurveySource4));
p4=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
    'DisplayName',[datestr(SurveyDatenum4) ' ' SurveySource4]);
% profile offset for location of dune design revetment and toe.
% revetment estimted to be from -1.5 to +6m from the peak elev location 
%  in the post construction truck survey. cobble toe +10.5m peak
[zmax,imax]=max(SM(SurvNum).Z1Dtransect);x0=SM(SurvNum).X1D(imax);

% 
% 
SurvNum=find([SM.Datenum] == SurveyDatenum6 & strcmp({SM.Source},SurveySource6));
p6=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
    'DisplayName',[datestr(SurveyDatenum6) ' ' SurveySource6]);

plot([-60 80],[.774 .774],'k--');text(-50,.65,'MSL','fontweight','bold');

p7=plot([x0+6 x0-1.5],[1.22 3.96],'k:','linewidth',5,'DisplayName','Buried Revetment 2T Stone');
p8=plot(x0+10.5,4,'m.','markersize',20,'DisplayName','Dune Design Toe');

set(gca,'xdir','reverse');grid on;
legend([p3 p1 p2 p4 p6 p7 p8],' 8 Apr 1998 Airborne LiDAR; post-El Nino eroded beach',...
    ' 6 Oct 2015  Airborne LiDAR; pre-El Nino BASELINE RECOVERED BEACH',...
    ' 3 Mar 2016  ATV-Dolly Survey; post-El Nino eroded beach',...
    '12 Sep 2019  Truck LiDAR; post-Dune Completion',...
    '9 Jan 2023  Truck LiDAR; most recent survey',...
    'Dune Design Buried Revetment 2T Stone','Dune Design Toe Location',...
    'location','northwest','fontsize',12);
set(gca,'xlim',[-60 80],'ylim',[0 9],'fontsize',12,'fontweight','demi')
xlabel('Mop Cross-Shore Profile Distance (m)')
ylabel('Elevation (m, NAVD88)')
title(['Mop ' num2str(MopNumber) ' Profiles'],'fontsize',16,'fontweight','bold')

print(gcf,'-dpng','-r100','-loose','DuneProfiles.png');
% 
% 
% 
% 
% 
% [X,Y,Z]=CombineSGgrids(669,679,SG(176).Datenum,SG(176).Source);
% figure;
% surf(X,Y,Z);shading flat;BeachColorbar;set(gca,'dataaspectratio',[1 1 .1]);view(2);
% for MopNumber=MopStart:MopEnd
% PlotLabelMopTransectUTM(MopNumber,'3d','k','ShadeOff');