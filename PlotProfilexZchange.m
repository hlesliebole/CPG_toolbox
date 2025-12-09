clear all
close all
addpath ..
MopNumber=582;
%nr=10; % number of most recent profiles to plot
nr=15; % number of most recent profiles to plot
% elevation levels Z to calculate dxZ/dt  
Z=0.5:.1:4.0;

%nr=2
figure('position',[66         192        1163         571]);
%  %----- add cobble sightings
%  
% matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
% load(matfile,'SA');
%  %  load Mop Transect Info
% load('MopTableUTM.mat','Mop');
% 
% % divide mop area into 20 mop subtransects at 1m xshore resolution,
% %  with an extra 100m of back beach for each
% [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);
% 
% idx=find(vertcat(SA.Class) > 1);
% Xutm=vertcat(SA.X);Yutm=vertcat(SA.Y);
% Xutm=Xutm(idx);Yutm=Yutm(idx);
% Z=vertcat(SA.Z);Z=Z(idx);
% dt=[];
% for n=1:size(SA,2)
%     dt=[dt' SA(n).Datenum*ones(size(SA(n).Z))']';
% end
% dt=dt(idx);
% 
% [dp,NearIdx]=...
%     pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);
% 
% [row,col] = ind2sub(size(xst),NearIdx);
% 
% hold on;pc=plot(x1d(col),Z,'m.','DisplayName','Past ATV Cobble Sightings');

%--------------

load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');

%figure;
nn=0;
xZ=[];
for n=1:nr    
    m=size(SM,2)-nr+n;
    dt(n)=SM(m).Datenum+.75; % add 18hrs to survey dates for ~6pm surveys
    if ~isnan(min(SM(m).Z1Dmean))
    nn=nn+1;
    z=SM(m).Z1Dmean;z(z <0)=NaN;
    
    nzl=0;
    for zl=Z
    nzl=nzl+1;
    xzl=intersections([SM(m).X1D(1) SM(m).X1D(end)],[zl zl],SM(m).X1D,SM(m).Z1Dmean);
    if isempty(xzl);xzl=NaN;end
    %fprintf('%6.1f %8.1f\n',zl,xzl)
    xZ(n,nzl)=xzl;
    end

%     p(nn)=plot(SM(m).X1D,z,'-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL,'%4.1f') 'm']);hold on;
    end
end

xZ=xZ';
dxZdt=diff(xZ,1,2);
dxZdtPcolor=horzcat(dxZdt,dxZdt(:,end));
dxZdtPcolor=vertcat(dxZdtPcolor,dxZdtPcolor(end,:));
%figure
pcolor(dt,min(Z)-0.05:.1:max(Z)+0.05,dxZdtPcolor)
set(gca,'color',[.5 .5 .5])
%pcolor(dxZdtPcolor)
polarmap;colormap(flipud(colormap))
set(gca,'fontsize',12)
cb=colorbar;
cb.Label.String='Contour Cross-shore Location Change Between Surveys (m)';
set(gca,'Xtick',dt,'xlim',[min(dt) max(dt)])
datetick('x','mm/dd HH:MM','keepticks');%set(gca,'xlim',[min(dt) max(dt)])
xtickangle(45)
xlabel('Local Time 2021');ylabel('Elevation (m, NAVD88) and Hs (m)')
title('Mop 582: Change in Xshore Locations of Subaerial Profile Contours Between Surveys')
set(gca,'fontsize',14)
set(cb,'fontsize',14)

hold on;
[tidetime,tidehgt]=getztide(floor(dt(1)),ceil(dt(end)),'lst','navd');
tidetime=tidetime-8*hours;
h_ax = gca;
h_ax_line = axes('position', get(h_ax, 'position'),'color','none');
pt=plot(tidetime,tidehgt,'k-','linewidth',2);
yl=get(h_ax,'ylim');
xl=get(h_ax,'xlim');
set(h_ax_line,'color','none','ylim',yl,'xtick',[],'ytick',[])
set(h_ax_line,'xlim',datetime(xl,'convertfrom','datenum'))


% hold on;
% 
h_ax2 = gca;
h_ax_line2 = axes('position', get(h_ax2, 'position'),'color','none');
[wavetime,wavehs]=getwaves(582,dt(1),dt(end),'lst');
%pw=plot(wavetime,wavehs,'b-');


ts1=timeseries(tidehgt,datenum(tidetime));
ts2=timeseries(wavehs,datenum(wavetime));
[tide,wave]=synchronize(ts1,ts2,'union');
pw=plot(wave,'b-','linewidth',1);hold on;
ptw=plot(tide+wave,'k:','linewidth',2);
datetick;
set(h_ax_line2,'color','none','ylim',yl,'xtick',[],'ytick',[],...
    'xlabel',[],'ylabel',[],'title',[]) 
set(h_ax_line2,'xlim',xl)
legend([pt pw ptw],'Tide','Hs','Tide+Hs','location','northwest')



% imagesc(flipud(dxZdt));shading flat;colorbar
% polarmap;colormap(flipud(colormap))
% % best historical fit index
% bfidx=279;
% nn=nn+1;
% z=SM(bfidx).Z1Dmean;z(z < -0.7)=NaN;
% p(nn)=plot(SM(bfidx).X1D,z,'k:','linewidth',2,'DisplayName',...
%         [datestr(SM(bfidx).Datenum,'mm/dd/yy') ' Best Past Fit z=0.5-1.0m']);
%     
% xl=get(gca,'xlim');
% set(gca,'xlim',[0 xl(2)]);
% %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
% 
% plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
% plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);
% ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
% set(gca,'xdir','reverse','fontsize',14);grid on;
% legend([p pc ps],'location','northwest');
% title(['MOP ' num2str(MopNumber) ' Mean Profiles']);
% xlabel('Cross-shore Distance (m)');
% ylabel('Elevation (m, NAVD88)');
% 
% makepng(['MOP' num2str(MopNumber) 'profiles.png'])
%     