
% uses intersections.m to find intersection of the profiles with different
%  elevation levels

close all
%clear all
%load M00582SM.mat
load M00553SM.mat
%SM([SM.Mopnum] ~= SM(1).Mopnum)=[];
% choose navd88 elevation that defines the shoreline location
ShoreElev=0.774; % MSL
%ShoreElev=1.344; % MHW
%ShoreElev=1.04; % halfway between MSL and MHW
% ShoreElev=1.566; % MHHW
% ShoreElev=2.119; % HAT
BarElev=-3.5;

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

% sort SM struct by Xshoreline location
  T=struct2table(SM); % sort by date before saving
  sortedT = sortrows(T, 'Xshoreline');
  %sortedT = sortrows(T, 'BeachVol');
  SM=table2struct(sortedT)';
  
  % find max nonNaN profile
  nsmax=find(~isnan([SM.Xshoreline]), 1, 'last' );

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
    if ~isnan(SM(ns).Xshoreline)
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
    if ~isnan(SM(ns).Xshoreline)
    plot(year(datetime(SM(ns).Datenum,'convert','datenum')),ns,'k.','markersize',10)
    end
end
axes(ax2);
set(ax2,'ylim',yl);set(gca,'xlim',[1 12]);
for ns=1:size(SM,2)
    if ~isnan(SM(ns).Xshoreline)
    plot(month(datetime(SM(ns).Datenum,'convert','datenum')),ns,'k.','markersize',10)
    end
end
axes(ax3);
set(ax3,'ylim',yl);set(ax3,'xlim',[0 .15],...
    'xtick',0:0.03:0.15,'xticklabels',[' 0 ';'.03';'.06';'.09';'.12';'.15']);
for ns=1:nsmax
    if ~isnan(SM(ns).Xshoreline)
    plot(SM(ns).ShorefaceSlope,ns,'k.','markersize',10)
    end
end
ax3t=axes('Position',ax3.Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold on;
set(ax3t,'ylim',yl);xlabel('Beach Width (m)','color','r');set(ax3t,'xcolor','r');
set(ax3t,'ylim',[1 nsmax],'ytick',[]);
for ns=1:nsmax
    if ~isnan(SM(ns).Xshoreline)
    plot(SM(ns).Xshoreline,ns,'ro','markersize',5)
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

%  Now reduce the grid to a 2d matrix of profiles 
%  with a 1m increment in the xshore location of the shoreline, or 
%  characteristic profiles and shoreface slope as a function 
%  of shoreline location.

xl=[0 xmax];
zl=[ShoreElev ShoreElev];
msl=[0.774 0.774];
mhhw=[1.566 1.566];
for i=1:size(Z,1)
       xshore(i)=max(intersections(xl,zl,0:xmax,Z(i,:)));
       xmsl(i)=max(intersections(xl,msl,0:xmax,Z(i,:)));
       xmhhw(i)=max(intersections(xl,mhhw,0:xmax,Z(i,:)));
       xslope(i)=(1.566-0.774)/abs(xmhhw(i)-xmsl(i));
end

% round shoreline locations to nearest 1m
xshore=round(xshore);
m=0;
C=nan(length(min(xshore):max(xshore)),xmax+1);
for n=min(xshore):max(xshore)
  m=m+1;
  CP(m).Xshoreline=n;
  idx=find(xshore == n); 
  CP(m).ShorefaceSlope=nanmean(xslope(idx));
  CP(m).X1D=0:xmax;
  CP(m).Z1Dmean=nanmean(Z(idx,:),1);
  C(m,:)=CP(m).Z1Dmean;
end

% make figure of data in charactersitic profile CP struct array

figure('position',[440 214 774 584]);

ax2=axes('position',[.25 .1 .5 .8]);
imagesc(0:xmax,min(xshore):max(xshore),C,'AlphaData',~isnan(C));
colormap(BeachColorMap);set(gca,'clim',[-8 5],'ydir','normal','xdir','reverse');
grid on;set(gca,'fontsize',13);
xlabel('Xshore Distance (m)');
title('Characterstic MOP Xshore Profiles vs. Mop Shoreline Location');

BeachBarColorbar

ax1=axes('position',[.06 .1 .15 .8]);
plot([CP.ShorefaceSlope],[CP.Xshoreline]);
ylabel('Mop MSL Shoreline X Location (m)');
xlabel({'MSL-MHHW';'Shoreface Slope'});
yl=get(ax2,'ylim');
set(gca,'ylim',yl,'fontsize',13);grid on;


