% example code plots all mediam profiles
%load M00552SM.mat
% close all
clearvars
Alat5=[];
Alon5=[];
Alat7=[];
Alon7=[];
%zoff=0;
zoff=-0.774;

for MopNumber=650
eval(['load M00' num2str(MopNumber) 'SM.mat'])
figure('position',[204 99 1120 637]);
hold on;



jumbo=find(contains({SM.File},'umbo'));
jumbo=find(contains({SM.File},'multibeam'));
clear p
for ns=1:size(SM,2)
    n=numel(SM(ns).X1D);
    if numel(SM(ns).Z1Dmedian) == n
        if ns == jumbo(end) %size(SM,2)
             %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dmedian+zoff,'m-','linewidth',2);
             p(1)=plot(SM(ns).X1D,SM(ns).Z1Dmedian+zoff,'b-','linewidth',3,'DisplayName',...
                 ['Latest Jumbo: ' datestr(SM(jumbo(end)).Datenum)]);
             n5m=find(SM(ns).Z1Dmedian+zoff < -5,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,SM(ns).X1D(n5m));
             % p(1)=plot(SM(ns).X1D(n5m),SM(ns).Z1Dmedian(n5m)+zoff,'ko','markersize',10,'markerfacecolor','g','DisplayName',...
             %     ['5m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat5=[Alat5 Lat];
             Alon5=[Alon5 Lon];
             n7m=find(SM(ns).Z1Dmedian+zoff < -7,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,SM(ns).X1D(n7m));
             % p(2)=plot(SM(ns).X1D(n7m),SM(ns).Z1Dmedian(n7m)+zoff,'ko','markersize',10,'markerfacecolor','y','DisplayName',...
             %     ['7m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat7=[Alat7 Lat];
             Alon7=[Alon7 Lon];

        elseif ns == jumbo(end)-1 %size(SM,2)
             %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dmedian+zoff,'m-','linewidth',2);
             p(2)=plot(SM(ns).X1D,SM(ns).Z1Dmedian+zoff,'-','color',[0 .75 0],'linewidth',3,'DisplayName',...
                 ['Previous Jumbo: ' datestr(SM(jumbo(end)-1).Datenum)]);

        elseif SM(ns).Datenum == 736440 %size(SM,2)
             %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dmedian+zoff,'m-','linewidth',2);
             p(3)=plot(SM(ns).X1D,SM(ns).Z1Dmedian+zoff,'r-','linewidth',3,'DisplayName',...
                 ['post-EN Jumbo: ' datestr(SM(ns).Datenum)]);hold on;

        else
    %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dmedian,'k-');
    p(4)=plot(SM(ns).X1D,SM(ns).Z1Dmedian+zoff,'-','Color',[.7 .7 .7],'DisplayName',...
                 ['All Profiles ']);hold on;
        end
    end
    hold on;
end
%view(3)
%datetick(gca,'y');
set(gca,'xdir','reverse','zlim',[-10 5]);
grid on;
xl=get(gca,'xlim');yl=get(gca,'ylim');
% fill3([xl xl(2) xl(1) xl(1)],[yl(1) yl yl(2) yl(1)],...
%     [.774 .774 .774 .774 .774],[.8 .8 .8],'FaceAlpha',0.9);
xlabel('Xshore Distance (m)');
%ylabel('Date');
ylabel('Elevation (m, MSL)');
title(['Mop ' num2str(SM(1).Mopnum) ' All Xshore Profiles']);
set(gca,'fontsize',14);
legend(p,'location','northwest','fontsize',18)
%makepng(['SolanaMop' num2str(MopNumber) 'Profiles.png'])
end

%save SolanaPUVsites.mat Alat5 Alon5 Alat7 Alon7