MopNumber=654;
load('M00654SA.mat')
[X1D,Z1D]=GetAllNearestPointsProfiles(SA,15,5);
X1D=X1D+15;
DateStart=datenum(2023,7,1);
DateEnd=datenum(2023,7,10);

ndx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx=ndx(end);

figure('position',[ 248   282   892   414]);
p1=plot(X1D,Z1D(ndx,:),'b-','linewidth',2);
hold on;


ndx=find((strcmp({SA.Source},'Gps') ) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx=ndx(end);
plot(X1D(137:end),Z1D(ndx,137:end),'b-','linewidth',2)



DateStart=datenum(2024,7,1);
DateEnd=datenum(2024,7,10);

ndx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx=ndx(end);
p2=plot(X1D,Z1D(ndx,:),'r-','linewidth',2);


ndx=find((strcmp({SA.Source},'Gps') ) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx=ndx(end);
plot(X1D(162:end),Z1D(ndx,162:end),'r-','linewidth',2)
xl=get(gca,'xlim');
%set(gca,'xlim',[0 xl(2)]);
set(gca,'xlim',[0 200]);xl(2)=200;
%yl=get(gca,'ylim');
set(gca,'ylim',[-4 6],'fontsize',16);

%plot(xl,[1.344 1.344],'k--');text(xl(2)-5,1.55,' MHW','fontsize',14,'fontweight','bold');
plot(xl,[1.566 1.566],'k--');text(xl(2)-5,1.766,' MHHW','fontsize',14,'fontweight','bold');
plot(xl,[.774 .774],'k:','linewidth',2);text(xl(2)-5,1.01,' MSL','fontsize',14,'fontweight','bold');
plot(xl,[-0.058 -0.058],'k--');text(xl(2)-5,0.16,' MLLW','fontsize',14,'fontweight','bold');
plot([48 48],[0.774 4.53],':','color','[.1 .1 .1]');
text(65,1.85,'\Leftarrow +40m','color','r','fontsize',16,'fontweight','bold');
text(72,1.05,'\Leftarrow +25m','color','r','fontsize',16,'fontweight','bold');
text(48,3.5,'\bf\Uparrow +3.7m (12ft)','color','r','fontsize',16,'fontweight','bold');
plot(68.5,1.566,'r.','markersize',20);plot(28,1.566,'b.','markersize',20)
plot(74,.774,'r.','markersize',20);plot(48,.774,'b.','markersize',20)
plot(48,4.53,'r.','markersize',20);

%ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
set(gca,'xdir','reverse','fontsize',16);grid on;
legend([p1 p2],'July 2023','July 2024','location','north','fontsize',16);
title(['Fletcher Cove Cross-Shore Profiles']);
xlabel('Cross-shore Distance (m)');
ylabel('Elevation (m, NAVD88)')
makepng('FletcherProfilesJuly2024.png')


