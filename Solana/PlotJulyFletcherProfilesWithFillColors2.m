MopNumber=654;
% load('M00654SA.mat')
% [X1D,Z1D]=GetAllNearestPointsProfiles(SA,15,5);
% X1D=X1D+15;

load FletcherProfileData.mat 

%% find and plot 2023 profile info

DateStart=datenum(2023,7,1);
DateEnd=datenum(2023,7,10);

ndx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx=ndx(end);

figure('position',[ 248   282   892   414]);
p1=plot(X1D,Z1D(ndx,:),'k-','linewidth',2);
hold on;

ndx1=find((strcmp({SA.Source},'Gps') ) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx1=ndx1(end);
plot(X1D(137:end),Z1D(ndx1,137:end),'k-','linewidth',2)

%% find and plot july 2024 profile info
DateStart=datenum(2024,7,1);
DateEnd=datenum(2024,7,10);

ndx2=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx2=ndx2(end);
p2=plot(X1D,Z1D(ndx2,:),'r-','linewidth',2);

%% fill in different profile areas with colors

ndx3=find((strcmp({SA.Source},'Gps') ) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);
ndx3=ndx3(end);

ixf=find(X1D > 80 & X1D < 169);
f1=fill([X1D(ixf) fliplr(X1D(ixf))],[Z1D(ndx1,ixf)  fliplr(Z1D(ndx3,ixf))],'y');%'facecolor',[0.9290 0.6940 0.1250]);
plot(X1D(137:end),Z1D(ndx1,137:end),'k-','linewidth',3)
plot(X1D(162:end),Z1D(ndx3,162:end),'r-','linewidth',3)

ixf=find(X1D >= 0 & X1D < 79);
xf=[X1D(ixf) fliplr(X1D(ixf))];
zf=[Z1D(ndx,ixf)  fliplr(Z1D(ndx2,ixf))];
igood=find(~isnan(zf));
f2=fill(xf(igood),zf(igood),'y');
set(f2,'facecolor',[0.9290 0.6940 0.1250])

ixf=find(X1D > 69 & X1D < 201);
xf=[X1D(ixf) 200 70];
zf=[max(vertcat(Z1D(ndx1,ixf),Z1D(ndx2,ixf),Z1D(ndx3,ixf))) 1.566 1.566] ;
f3=fill(xf,zf,'c');%'facecolor',[0.9290 0.6940 0.1250]);

%% replot the profile lines

plot(X1D(137:end),Z1D(ndx,137:end),'k-','linewidth',3)
plot(X1D(137:end),Z1D(ndx1,137:end),'k-','linewidth',3)
plot(X1D(162:end),Z1D(ndx2,162:end),'r-','linewidth',3)
plot(X1D(162:end),Z1D(ndx3,162:end),'r-','linewidth',3)

p1=plot(X1D,Z1D(ndx,:),'k-','linewidth',3);
p2=plot(X1D,Z1D(ndx2,:),'r-','linewidth',3);

%% add water level lines
xl=get(gca,'xlim');
set(gca,'xlim',[0 200]);xl(2)=200;
set(gca,'ylim',[-4 6],'fontsize',16);

plot(xl,[1.566 1.566],'k--');text(xl(2)-5,1.766,' MHHW','fontsize',14,'fontweight','bold');
plot(xl,[.774 .774],'k:','linewidth',2);text(xl(2)-5,1.01,' MSL','fontsize',14,'fontweight','bold');
plot(xl,[-0.058 -0.058],'k--');text(xl(2)-5,0.16,' MLLW','fontsize',14,'fontweight','bold');

%% add arrows and labels showing change
plot([48 48],[0.774 4.53],':','color','[.1 .1 .1]');
text(65,1.85,'\Leftarrow +40m','color','k','fontsize',16,'fontweight','bold');
text(72,1.05,'\Leftarrow +25m','color','k','fontsize',16,'fontweight','bold');
text(48,3.5,'\bf\Uparrow +3.7m (12ft)','color','k','fontsize',16,'fontweight','bold');
plot(68.5,1.566,'r.','markersize',20);plot(28,1.566,'b.','markersize',20)
plot(74,.774,'r.','markersize',20);plot(48,.774,'b.','markersize',20)
plot(48,4.53,'r.','markersize',20);

%% add legend
legend([p1 p2 f2 f1],'July 2023','July 2024','Coarse Nourishment Sand',...
    'Coarse + Natural Fine Sand Mix','location','northwest','fontsize',16);

%% reverse and lable axis
set(gca,'xdir','reverse','fontsize',16);grid on;
title('Fletcher Cove Cross-Shore Profiles','fontsize',20);
xlabel('Cross-shore Distance (m)');
ylabel('Elevation (m, NAVD88)')

%% make png image of figure
set(gcf,'Units','pixels','inverthardcopy','off');
scrpos = get(gcf,'Position');
newpos = scrpos/90;
set(gcf,'PaperUnits','inches',...
    'PaperPosition',newpos)
print(gcf,'-dpng','-r90','FletcherProfilesJuly2024.png')

