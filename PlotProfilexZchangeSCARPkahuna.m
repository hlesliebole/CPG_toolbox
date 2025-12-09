clear all
close all
addpath ..
MopNumber=582;
%nr=10; % number of most recent profiles to plot
nr=15; % number of most recent profiles to plot
% elevation levels Z to calculate dxZ/dt  
Z=-.2:.1:4.0;

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

S1=datenum(2000,10,1);
S2=datenum(2001,9,30);
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');


idx=find([SM.Datenum] >= S1 & [SM.Datenum] <= S2);
[ud,idu]=unique([SM(idx).Datenum]);
idx=idx(idu); % one survey per day

%figure;
nn=0;
xZ=[];
nr=numel(idx);
for n=1:nr    
    m=idx(n);
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
    xZ(n,nzl)=xzl(1);
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
pcolor(dt,min(Z)-0.05:.1:max(Z)+0.05,dxZdtPcolor);
set(gca,'clim',[-8 8]);
set(gca,'color',[.8 .8 .8])
%pcolor(dxZdtPcolor)
polarmap;colormap(flipud(colormap))
set(gca,'fontsize',12)
cb=colorbar;
cb.Label.String='Contour Cross-shore Location Change Between Surveys (m)';
set(gca,'Xtick',dt,'xlim',[min(dt) max(dt)])
datetick('x','mm/dd/yy');set(gca,'xlim',[min(dt) max(dt)])
xtickangle(30)
xlabel('Local Time');ylabel({'Elevation (m, NAVD88) and Hs (m)',' '})
title({[datestr(dt(1),'mm/dd/yyyy') ' to ' datestr(dt(end),'mm/dd/yyyy')],...
    'Mop 582: Change in Xshore Locations of Subaerial Profile Contours Between Surveys'})
set(gca,'fontsize',14)
set(cb,'fontsize',14)

hold on;
[tidetime,tidehgt]=getztide(floor(dt(1)),ceil(dt(end)),'lst','navd');
tidetime=tidetime-8*hours;
% 
h_ax = gca;
h_ax_line = axes('position', get(h_ax, 'position'),'color','none');
% 
% pt=plot(tidetime,tidehgt,'k-','linewidth',2);
hold on
% add tide levels
plot(tidetime,tidehgt*0-0.058,'g-','linewidth',2);
ptime=datetime(dt(1),'convertfrom','datenum');
text(ptime,-0.1,'MLLW ','horizontalalign','right',...
    'fontsize',12,'fontweight','bold','color',[0 .8 0])
plot(tidetime,tidehgt*0+0.774,'g-','linewidth',2);
text(ptime,.774,'MSL ','horizontalalign','right',...
    'fontsize',12,'fontweight','bold','color',[0 .8 0])
plot(tidetime,tidehgt*0+1.566,'g-','linewidth',2);
text(ptime,1.65,'MHHW ','horizontalalign','right',...
    'fontsize',12,'fontweight','bold','color',[0 .8 0])
% 
yl=get(h_ax,'ylim');
xl=get(h_ax,'xlim');
set(h_ax_line,'color','none','ylim',yl,'xtick',[],'ytick',[])
set(h_ax_line,'xlim',datetime(xl,'convertfrom','datenum'))
% 
% 
% % hold on;
% % 
h_ax2 = gca;
h_ax_line2 = axes('position', get(h_ax2, 'position'),'color','none');
[wavetime,wavehs]=getwaves(582,dt(1),dt(end),'lst');
% %pw=plot(wavetime,wavehs,'b-');
% 
% 
ts1=timeseries(tidehgt,datenum(tidetime));
ts2=timeseries(wavehs,datenum(wavetime));
[tide,wave]=synchronize(ts1,ts2,'union');

%plot(wave,'-','color',[.4 .4 .4],'linewidth',3);hold on;
%pw=plot(wave,'k-','linewidth',1);hold on;
%ptw=plot(tide+wave,'k:','linewidth',1);

maxhs=wavehs;
maxtide=tidehgt;
mintide=tidehgt;
yd=day(tidetime,'dayofyear');
%ydw=day(wave.Time,'dayofyear');
ydtw=day(datetime([wave.Time],'convertfrom','datenum'),'dayofyear');
maxtwave=tide+wave;
for iyd=1:365
    idx=find(yd == iyd);
    if ~isempty(idx)
        maxtide(idx)=max(tidehgt(idx));
        mintide(idx)=min(tidehgt(idx));
    end
    idx=find(ydtw == iyd);
    if ~isempty(idx)
        maxtwave.Data(idx)=max(maxtwave.Data(idx));
    end
end
ts3=timeseries(maxtide,datenum(tidetime));
ts4=timeseries(mintide,datenum(tidetime));
ts5=timeseries(maxhs,datenum(wavetime));
[tide,wave]=synchronize(ts3,ts5,'union');
plot(ts3,'-','color',[.4 .4 .4],'linewidth',3);hold on;
pmt=plot(ts3,'g-','linewidth',1);hold on;
plot(ts4,'-','color',[.4 .4 .4],'linewidth',3);hold on;
pmt2=plot(ts4,'y-','linewidth',1);hold on;
plot(maxtwave,'-','color',[.4 .4 .4],'linewidth',3);hold on;
pw=plot(maxtwave,'m-','linewidth',1);hold on;


datetick;
set(h_ax_line2,'xlim',xl,'ylim',yl,'color','none','xtick',[],'ytick',[],...
    'xlabel',[],'ylabel',[],'title',[]) 
set(h_ax_line2,'xlim',xl)
legend([pw pmt pmt2],'Max Daily Tide+Hs','Max Daily Tide','Min Daily Tide','location','northwest')
set(gcf,'inverthardcopy','off')
% 
% 
% 
% % imagesc(flipud(dxZdt));shading flat;colorbar
% % polarmap;colormap(flipud(colormap))
% % % best historical fit index
% % bfidx=279;
% % nn=nn+1;
% % z=SM(bfidx).Z1Dmean;z(z < -0.7)=NaN;
% % p(nn)=plot(SM(bfidx).X1D,z,'k:','linewidth',2,'DisplayName',...
% %         [datestr(SM(bfidx).Datenum,'mm/dd/yy') ' Best Past Fit z=0.5-1.0m']);
% %     
% % xl=get(gca,'xlim');
% % set(gca,'xlim',[0 xl(2)]);
% % %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
% % 
% % plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
% % plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);
% % ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
% % set(gca,'xdir','reverse','fontsize',14);grid on;
% % legend([p pc ps],'location','northwest');
% % title(['MOP ' num2str(MopNumber) ' Mean Profiles']);
% % xlabel('Cross-shore Distance (m)');
% % ylabel('Elevation (m, NAVD88)');
% % 
% % makepng(['MOP' num2str(MopNumber) 'profiles.png'])
% %     