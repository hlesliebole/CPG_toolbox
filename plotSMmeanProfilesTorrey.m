% example code plots all mediam profiles
%load M00552SM.mat
close all
clearvars
Alat5=[];
Alon5=[];
Alat7=[];
Alon7=[];
Alat10=[];
Alon10=[];
Alat15=[];
Alon15=[];
%zoff=0;
zoff=-0.774;


for MopNumber=[580 586]
%eval(['load M00' num2str(MopNumber) 'SM.mat'])
eval(['load M00' num2str(MopNumber) 'SA.mat'])
[X1D,Z1D]=GetNearestPointsProfiles(SA,25,5);

figure('position',[204 99 1120 637]);
hold on;

jumbo=find(contains({SA.File},'umbo'));
clear p
idone=0;
for ns=1:size(Z1D,1)
    n=numel(X1D);
    if numel(Z1D(ns,:)) == n
        if ns == jumbo(end) %size(SM,2)
             %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dtransect+zoff,'m-','linewidth',2);
             p(5)=plot(X1D,Z1D(ns,:)+zoff,'b-','linewidth',2,'DisplayName',...
                 ['Latest Jumbo: ' datestr(SA(jumbo(end)).Datenum)]);
             n5m=find(Z1D(ns,:)+zoff < -5,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X1D(n5m));
             p(1)=plot(X1D(n5m),Z1D(ns,n5m)+zoff,'ko','markersize',10,'markerfacecolor','g','DisplayName',...
                 ['5m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat5=[Alat5 Lat];
             Alon5=[Alon5 Lon];
             n7m=find(Z1D(ns,:)+zoff < -7,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X1D(n7m));
             p(2)=plot(X1D(n7m),Z1D(ns,n7m)+zoff,'ko','markersize',10,'markerfacecolor','y','DisplayName',...
                 ['7m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat7=[Alat7 Lat];
             Alon7=[Alon7 Lon];
             n10m=find(Z1D(ns,:)+zoff < -10,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X1D(n10m));
             p(3)=plot(X1D(n10m),Z1D(ns,n10m)+zoff,'ko','markersize',10,'markerfacecolor','c','DisplayName',...
                 ['10m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat10=[Alat10 Lat];
             Alon10=[Alon10 Lon];


        % elseif SM(ns).Datenum == 736440 %size(SM,2)
        %      %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dtransect+zoff,'m-','linewidth',2);
        %      p(5)=plot(SM(ns).X1D,SM(ns).Z1Dtransect+zoff,'r-','linewidth',2,'DisplayName',...
        %          ['post-EN Jumbo: ' datestr(SM(ns).Datenum)]);hold on;

        else
    %plot3(SM(ns).X1D,SM(ns).Datenum*ones(1,n),SM(ns).Z1Dtransect,'k-');
      plot(X1D,Z1D(ns,:)+zoff,'-','Color',[.7 .7 .7]);
             
             %if min(Z1D(ns,:)) < -15 & idone == 0
              if min(Z1D(ns,:)) < -15 & strcmp(SA(ns).Source,'USACE') & idone == 0
              n15m=find(Z1D(ns,:)+zoff < -15,1,'first');
             [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X1D(n15m));
             p(4)=plot(X1D(n15m),Z1D(ns,n15m)+zoff,'ko','markersize',10,'markerfacecolor','r','DisplayName',...
                 ['15m point Lat:' num2str(Lat,'%9.5f') '  Lon: ' num2str(Lon,'%9.5f')]);
             Alat15=[Alat15 Lat];
             Alon15=[Alon15 Lon];
             idone = 1;
             end

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
title(['Mop ' num2str(SA(1).Mopnum) ' All Xshore Profiles']);
set(gca,'fontsize',14);
legend(p,'location','northwest','fontsize',18)
makepng(['TorreyMop' num2str(MopNumber) 'Profiles.png'])
end

save TorreyPUVsites.mat Alat5 Alon5 Alat7 Alon7 Alat10 Alon10 Alat15 Alon15