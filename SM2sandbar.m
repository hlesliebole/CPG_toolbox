
% uses intersections.m to find intersection of the profiles with different
%  elevation levels

close all
clear all
% MopNumber=582;
% load M00582SM.mat
MopNumber=553;
load M00553SM.mat
%SM([SM.Mopnum] ~= SM(1).Mopnum)=[];
% choose navd88 elevation that defines the shoreline location
ShoreElev=0.774; % MSL
%ShoreElev=1.344; % MHW
%ShoreElev=1.04; % halfway between MSL and MHW
% ShoreElev=1.566; % MHHW
% ShoreElev=2.119; % HAT
BarElev=-5.0;

% loop through mean profiles and add the shoreline x location
%  to the struct array
   
for ns=1:size(SM,2)
    % remove any negative xshore profile values
    SM(ns).Z1Dmean(SM(ns).X1D < 0)=NaN;
    % find shoreline datum elevation
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)];
    zl=[ShoreElev ShoreElev];   
    xz=intersections(xl,zl,SM(ns).X1D,SM(ns).Z1Dmean);
    if isempty(xz)
        SM(ns).Xshoreline=NaN;
        SM(ns).BeachVol=NaN;
    else
        subair=find(SM(ns).X1D > 0 & SM(ns).X1D <= max(xz));
        SM(ns).Xshoreline=max(xz);
        ixd=find(~isnan(SM(ns).Z1Dmean));
        zi=interp1(SM(ns).X1D(ixd),SM(ns).Z1Dmean(ixd),SM(ns).X1D,'linear','extrap');
        SM(ns).BeachVol=nansum(zi(subair));
    end
    
    % find sandbar datum elevation
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)];
    zl=[BarElev BarElev];   
    xz=intersections(xl,zl,SM(ns).X1D,SM(ns).Z1Dmean);
    if isempty(xz)
        SM(ns).Xbar=NaN;
    else
        SM(ns).Xbar=max(xz);
    end
    
    % find MSL-MHW shoreface slope
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)];
     
    zl=[0.774 0.774];
    x1=intersections(xl,zl,SM(ns).X1D,SM(ns).Z1Dmean);
    if ~isempty(x1); x1=nanmin(x1); end
    zl=[1.566 1.566];
    x2=intersections(xl,zl,SM(ns).X1D,SM(ns).Z1Dmean);
    if ~isempty(x2); x2=nanmin(x2); end
    
    if ~isempty(x1) && ~isempty(x2)
        SM(ns).ShorefaceSlope=(1.566-0.774)/abs(x2-x1);
    else
        SM(ns).ShorefaceSlope=NaN;
        % also eliminate shoreline location and beach vol if insufficient data
        %  for the shoreface slope calc
        SM(ns).Xshoreline=NaN;
        SM(ns).BeachVol=NaN;
    end

end

% sort SM struct by Xbar location
  T=struct2table(SM); % sort by date before saving
  sortedT = sortrows(T, 'Xbar');
  %sortedT = sortrows(T, 'BeachVol');
  SM=table2struct(sortedT)';
  
  % find max nonNaN profile
  nsmax=find(~isnan([SM.Xbar]), 1, 'last' );

figure('position',[204 99 1120 637]);

ax1=axes('position',[.05 .1 .12 .8],'fontsize',14);hold on;
grid on;box on;xlabel('Year');set(gca,'ylim',[1 nsmax]);
ylabel('Shoreline Location Rank Number');
ax2=axes('position',[.18 .1 .12 .8],'fontsize',14);hold on;
grid on;box on;set(gca,'ylim',[1 nsmax]);
set(gca,'yticklabels',[]');xlabel('Month');
ax3=axes('position',[.32 .1 .12 .8],'fontsize',14);hold on;
grid on;box on;set(gca,'ylim',[1 nsmax]);
set(gca,'yticklabels',[]');xlabel({'MSL-MHHW';'Shoreface Slope'});

ax4=axes('position',[.5 .1 .4 .8]);hold on;set(gca,'ylim',[1 nsmax]);

load BeachBarColorMap
cm = BeachColorMap;
    
for ns=1:nsmax
    if ~isnan(SM(ns).Xbar)
    n=length(SM(ns).X1D);
    x=SM(ns).X1D;
    %y=SM(ns).Datenum*ones(1,n);
    %y=SM(ns).Xshoreline*ones(1,n);
    y=ns*ones(1,n);
    z=SM(ns).Z1Dmean;
    x(isnan(z))=[];
    y(isnan(z))=[];
    z(isnan(z))=[];   
    
    zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    (BeachColorBarLims(2)-BeachColorBarLims(1));
    zscaled(zscaled < 1)=1;
    zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);   

    scp=scatter3(x, y, z, 10, cm(ceil(zscaled),:), 'filled');
    scp.MarkerFaceAlpha = .9;
    scp.MarkerEdgeAlpha = .9;
     
    hold on;
    zl=[BarElev BarElev];
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)];
    xz=intersections(xl,zl,SM(ns).X1D,SM(ns).Z1Dmean);
    %plot3(xz,ns*ones(size(xz)),BarElev*ones(size(xz))+.25,'k.','markersize',15);
    
    end  
end

view(2) % top view

%datetick(gca,'y');
set(gca,'xdir','reverse')
%set(gca,'ydir','reverse')
grid on;

xl=get(gca,'xlim');yl=get(gca,'ylim');
% fill3([xl xl(2) xl(1) xl(1)],[yl(1) yl yl(2) yl(1)],...
%     [.774 .774 .774 .774 .774],[.8 .8 .8],'FaceAlpha',0.5);
xlabel('Xshore Distance (m)');
%ylabel('Ranked Shoreline Location');
zlabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(SM(1).Mopnum) ' Area Mean Xshore Profiles']);
set(gca,'fontsize',14);

BeachBarColorbar

yl=get(ax4,'ylim');
set(ax1,'ylim',yl);
axes(ax1);
for ns=1:nsmax
    if ~isnan(SM(ns).Xbar)
    plot(year(datetime(SM(ns).Datenum,'convert','datenum')),ns,'k.','markersize',10)
    end
end
axes(ax2);
set(ax2,'ylim',yl);set(gca,'xlim',[1 12]);
for ns=1:size(SM,2)
    if ~isnan(SM(ns).Xbar)
    plot(month(datetime(SM(ns).Datenum,'convert','datenum')),ns,'k.','markersize',10)
    end
end
axes(ax3);
set(ax3,'ylim',yl);set(ax3,'xlim',[0 .15],...
    'xtick',0:0.03:0.15,'xticklabels',[' 0 ';'.03';'.06';'.09';'.12';'.15']);
for ns=1:nsmax
    if ~isnan(SM(ns).Xbar)
    plot(SM(ns).ShorefaceSlope,ns,'k.','markersize',10)
    end
end
ax3t=axes('Position',ax3.Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold on;
set(ax3t,'ylim',yl);xlabel('Beach Width (m)','color','r');set(ax3t,'xcolor','r');
set(ax3t,'ylim',[1 nsmax],'ytick',[]);
for ns=1:nsmax
    if ~isnan(SM(ns).Xbar)
    plot(SM(ns).Xbar,ns,'ro','markersize',5)
    end
end

%  make 2d grid of ranked profiles
xmax=max([SM(1:nsmax).X1D]);  % max crossshore location
[X,Y]=meshgrid(0:xmax,1:nsmax);
N=[];Xg=[];Zg=[];
for ns=1:nsmax;X1=SM(ns).X1D;Z1=SM(ns).Z1Dmean;...
        X1(isnan(Z1))=[];Z1(isnan(Z1))=[];...
        Xg=[Xg X1];Zg=[Zg Z1];N=[N X1*0+ns];
end
Z=griddata(Xg,N,Zg,X,Y);
Z=movmean(Z,round(nsmax/6),1); % smooth on ~2 month shoreline position time scale
figure;surf(X,Y,Z);shading flat;set(gca,'xdir','reverse');
colormap(BeachColorMap);set(gca,'clim',[-8 5]);view(2)

figure;histogram([SM.Xbar]);xlabel('Barface Xshore Location (m)');ylabel('N obs');

nf=0;
idx=find(~isnan([SM.Xbar]) & month(datetime([SM.Datenum],'convertfrom','datenum')) > 3 & ...
month(datetime([SM.Datenum],'convertfrom','datenum')) < 10);idx=fliplr(idx);
figure('position',[201 79 1000 718]);hold on;
p1=plot(SM(idx(1)).X1D,SM(idx(1)).Z1Dmean,'b-','linewidth',2);
p2=plot(SM(idx(1)).X1D,SM(idx(1)).Z1Dmean,'b-','linewidth',2);
p3=plot(SM(idx(1)).X1D,SM(idx(1)).Z1Dmean,'b-','linewidth',2);
set(gca,'xdir','reverse','ylim',[-8 5]);
xl=get(gca,'xlim');set(gca,'xlim',xl);
title([{['Mop: ' num2str(MopNumber)]},...
    {[' Observed April to September Profile Evolution Relative to ' num2str(BarElev) 'm Contour Xshore Location']}],...
    'fontsize',16);
t1=text(mean(xl),3,datestr(SM(idx(1)).Datenum),'horizontalalign','center',...
    'fontweight','bold','fontsize',20);
xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
grid on;box on;set(gca,'fontsize',12)
plot(xl,[BarElev BarElev],'k:','linewidth',2);
plot([min([SM(idx).Xbar]) max([SM(idx).Xbar])],[BarElev BarElev],'k-','linewidth',3);

plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',13);
plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',13);
plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',13);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',13);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',13);
nf=nf+1;M(nf)=getframe(gcf);
for n=2:numel(idx)
    %pause
    delete(p1);delete(p2);delete(p3);
  if n > 2
    p1=plot(SM(idx(n-2)).X1D,SM(idx(n-2)).Z1Dmean,'-','linewidth',2,'color',[.9 .9 .9]);
    p2=plot(SM(idx(n-1)).X1D,SM(idx(n-1)).Z1Dmean,'-','linewidth',2,'color',[.7 .7 .7]);
    p3=plot(SM(idx(n)).X1D,SM(idx(n)).Z1Dmean,'b-','linewidth',2);
  else      
    p1=plot(SM(idx(n-1)).X1D,SM(idx(n-1)).Z1Dmean,'-','linewidth',2,'color',[.7 .7 .7]);
    p2=plot(SM(idx(n-1)).X1D,SM(idx(n-1)).Z1Dmean,'-','linewidth',2,'color',[.7 .7 .7]);
    p3=plot(SM(idx(n)).X1D,SM(idx(n)).Z1Dmean,'b-','linewidth',2);   
  end
  delete(t1)
  t1=text(mean(xl),3,datestr(SM(idx(n)).Datenum),'horizontalalign','center',...
    'fontweight','bold','fontsize',20);
   nf=nf+1;M(nf)=getframe(gcf);
end

v=VideoWriter('Mop553ProfileEvolution.mp4','MPEG-4');
v.FrameRate=1;
open(v)
writeVideo(v,M)
close(v)

% for xbar=round(min([SM.Xbar])):5:round(max([SM.Xbar]))
%  [xmin,imin]=min(abs([SM(idx).Xbar]-xbar));
%  plot(SM(idx(imin)).X1D,SM(idx(imin)).Z1Dmean,'-');
% end



%  [xmin,imin]=min(abs([SM(idx).Xbar]-xbar));
%  plot(SM(idx(imin)).X1D,SM(idx(imin)).Z1Dmean,'-');
%  plot(SM(idx(end)).X1D,SM(idx(end)).Z1Dmean,'-');set(gca,'xdir','reverse');

% %  Now reduce the grid to a 2d matrix of profiles 
% %  with a 1m increment in the xshore location of the shoreline, or 
% %  characteristic profiles and shoreface slope as a function 
% %  of shoreline location.
% 
% xl=[0 xmax];
% zl=[ShoreElev ShoreElev];
% msl=[0.774 0.774];
% mhhw=[1.566 1.566];
% for i=1:size(Z,1)
%        xshore(i)=max(intersections(xl,zl,0:xmax,Z(i,:)));
%        xmsl(i)=max(intersections(xl,msl,0:xmax,Z(i,:)));
%        xmhhw(i)=max(intersections(xl,mhhw,0:xmax,Z(i,:)));
%        xslope(i)=(1.566-0.774)/abs(xmhhw(i)-xmsl(i));
% end
% 
% % round shoreline locations to nearest 1m
% xshore=round(xshore);
% m=0;
% C=nan(length(min(xshore):max(xshore)),xmax+1);
% for n=min(xshore):max(xshore)
%   m=m+1;
%   CP(m).Xshoreline=n;
%   idx=find(xshore == n); 
%   CP(m).ShorefaceSlope=nanmean(xslope(idx));
%   CP(m).X1D=0:xmax;
%   CP(m).Z1Dmean=nanmean(Z(idx,:),1);
%   C(m,:)=CP(m).Z1Dmean;
% end
% 
% % make figure of data in charactersitic profile CP struct array
% 
% figure('position',[440 214 774 584]);
% 
% ax2=axes('position',[.25 .1 .5 .8]);
% imagesc(0:xmax,min(xshore):max(xshore),C,'AlphaData',~isnan(C));
% colormap(BeachColorMap);set(gca,'clim',[-8 5],'ydir','normal','xdir','reverse');
% grid on;set(gca,'fontsize',13);
% xlabel('Xshore Distance (m)');
% title('Characterstic MOP Xshore Profiles vs. Mop Shoreline Location');
% 
% BeachBarColorbar
% 
% ax1=axes('position',[.06 .1 .15 .8]);
% plot([CP.ShorefaceSlope],[CP.Xshoreline]);
% ylabel('Mop MSL Shoreline X Location (m)');
% xlabel({'MSL-MHHW';'Shoreface Slope'});
% yl=get(ax2,'ylim');
% set(gca,'ylim',yl,'fontsize',13);grid on;
% 
% 
