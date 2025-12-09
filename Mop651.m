load M00651SM.mat

idx=find(strcmp({SM.Source},'Multibeam'));
ndx=find(contains({SM.File},'umbo'));
ndx=ndx(end);

figure('position',[ 126         114        1069         630]);
subplot(2,1,1)
plot(SM(1).X1D,SM(idx(1)).Z1Dtransect,'k-','displayname',...
    datestr(SM(idx(1)).Datenum));
set(gca,'xdir','reverse','fontsize',16)
title('Mop 651 Transect Line Multibeam Data')
grid on
legend('location','northwest')
hold on;
plot(SM(1).X1D,SM(idx(2)).Z1Dtransect,'r-','displayname',...
    datestr(SM(idx(2)).Datenum));
plot(SM(1).X1D,SM(ndx).Z1Dtransect,'b-','displayname',...
    [ datestr(SM(ndx).Datenum) ' Jumbo']);
xlabel('Xshore Distance (m)')
ylabel('Elevation (m, NAVD88)')
set(gca,'xlim',[0 900])
subplot(2,1,2)
plot(SM(1).X1D,SM(idx(2)).Z1Dtransect-SM(idx(1)).Z1Dtransect,'r-');
hold on;
plot(SM(1).X1D,SM(ndx).Z1Dtransect-SM(idx(1)).Z1Dtransect,'b-');
hold on;
set(gca,'xdir','reverse','fontsize',16)
grid on
xlabel('Xshore Distance (m)')
title('Difference from 21 Nov Multibeam')
ylabel('Difference (m)')
set(gca,'xlim',[0 900])
makepng('Mop651MultibeamDiffNovDec.png')