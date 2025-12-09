
% ContourTracker.m

% Plots possible target contour migration from its location in the most recent
%  survey in the CPG MOP database, based on historical profiles.

% uses intersections.m to find intersection of the profiles with different
%  elevation levels
clear all
% settings
MopNumber=582; % Mop number
%Tcontour=-0.058;% MLLW % navd88 depth of contour to try and track
%ShoreElev=0.774; % MSL % navd88 elevation that defines the shoreline location
ShoreElev=.774; 
dx=7; % +/- shoreline location change plotting tolerance (m) 

close all

figure('Position',[190   192   994   532]);
 %----- add cobble sightings
 
matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');
 %  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);

idx=find(vertcat(SA.Class) > 1);
Xutm=vertcat(SA.X);Yutm=vertcat(SA.Y);
Xutm=Xutm(idx);Yutm=Yutm(idx);
Z=vertcat(SA.Z);Z=Z(idx);
dt=[];
for n=1:size(SA,2)
    dt=[dt' SA(n).Datenum*ones(size(SA(n).Z))']';
end
dt=dt(idx);

[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

[row,col] = ind2sub(size(xst),NearIdx);

hold on;pc=plot(x1d(col),Z,'m.');

%--------------


% load survey morphology cpg mop file
SMmatfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
load(SMmatfile,'SM');

% loop through mean profiles and add the shoreline x location,
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
        %SM(ns).BeachVol=NaN;
    else
        subair=find(SM(ns).X1D > 0 & SM(ns).X1D <= max(xz));
        SM(ns).Xshoreline=max(xz);
%         ixd=find(~isnan(SM(ns).Z1Dmean));
%         zi=interp1(SM(ns).X1D(ixd),SM(ns).Z1Dmean(ixd),SM(ns).X1D,'linear','extrap');
%         SM(ns).BeachVol=nansum(zi(subair));
    end
   
end

% sort SM struct by Xshoreline location in new struct array SSM
  T=struct2table(SM); % sort by date before saving
  sortedT = sortrows(T, 'Xshoreline');
  SSM=table2struct(sortedT)'; % narrowest-widest shoreline struct array
  
  
  % most recent shoreline x location
  mrX=SM(end).Xshoreline;
  
%   % index of surveys with same shoreline location +/- dx setting
%   idx=find(vertcat(SM.Xshoreline) >= mrX-dx & vertcat(SM.Xshoreline) <= mrX+dx);
  
  % index of surveys with same shoreline location +/- dx setting
  idx=find([SM.Datenum] >= datenum(2019,10,13) & [SM.Datenum] <= datenum(2020,4,1));
  
  clear M
  for nn=1:numel(idx)
      if nn == 1
          p(nn)=plot(SM(idx(nn)).X1D,SM(idx(nn)).Z1Dmean,'k-','linewidth',2);hold on;
  
  set(gca,'ylim',[-1 2.5]);
  xl=get(gca,'xlim');if xl(1) < 0; xl(1) = 0; set(gca,'xlim',[0 120]);end
  set(gca,'fontsize',14);
  set(gca,'xtick',[0:5:500]);
  set(gca,'xdir','reverse','xgrid','on');

  % lable misc tide elevations a p sensor elev
  ztide=[-0.631 -.058 .218 .774 1.344 1.566 2.119 ];
  for n=1:length(ztide)
              
                 plot(xl,[ztide(n) ztide(n)],'k:','linewidth',2);
                 plot(xl,[ztide(n) ztide(n)],'b:','linewidth',1); 
  end
        plot(xl,[ztide(2)-.5 ztide(2)-.5],'r--');
        text(xl(end),ztide(1),' LAT','fontsize',14,'verticalalign','bottom');
        text(79,ztide(2)-0.5,' P','color','r','verticalalign','bottom','fontweight','bold');
        text(xl(end),ztide(2),' MLLW','verticalalign','bottom');
        text(xl(end),ztide(3),' MLW','verticalalign','bottom');
        text(xl(end),ztide(4),' MSL','fontsize',14,'verticalalign','bottom');
        text(xl(end),ztide(5),' MHW','verticalalign','bottom');
        text(xl(end),ztide(6),' MHHW','verticalalign','bottom');
        text(xl(end),ztide(7),' HAT','fontsize',14,'verticalalign','bottom');
        
        %plot(SM(end).Xshoreline,ShoreElev,'k.','markersize',20)
        
%  legend([pl pe],datestr(SM(end).Datenum,'mm/dd/yyyy'),'Most Eroded',...
%      'location','northwest')
 
 xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88');
 
 title(['Historical ' num2str(MopNumber) ' Profile Evolution from 14 Oct 2019']);
          
          tx=text(80,1.8,datestr(SM(idx(nn)).Datenum),'fontsize',35);
          M(nn)=getframe(gcf);
      else
          if nn > 2;set(p(nn-1),'color',[.8 .8 .8],'linewidth',1);end
          p(nn)=plot(SM(idx(nn)).X1D,SM(idx(nn)).Z1Dmean,'-','linewidth',2);hold on;
          set(gca,'xlim',[0 120])
          delete(tx)
          tx=text(80,1.8,datestr(SM(idx(nn)).Datenum),'fontsize',35);
          text(79,ztide(2)-0.5,' P','color','r','verticalalign','bottom','fontweight','bold');
          M(nn)=getframe(gcf);
      end
  end
  
v=VideoWriter('ProfileEvolution.mp4','MPEG-4');
v.FrameRate=1;
open(v)
writeVideo(v,M)
close(v)
%   
%   pl=plot(SM(end).X1D,SM(end).Z1Dmean,'r-','linewidth',2);hold on;
%   %pe=plot(SSM(1).X1D,SSM(1).Z1Dmean,'k--','linewidth',2);hold on;
%   set(gca,'ylim',[-1 2.5]);
%   xl=get(gca,'xlim');if xl(1) < 0; xl(1) = 0; set(gca,'xlim',xl);end
%   set(gca,'fontsize',14);
%   set(gca,'xtick',[0:5:500]);
%   set(gca,'xdir','reverse','xgrid','on');
% 
%   % lable misc tide elevations a p sensor elev
%   ztide=[-0.631 -.058 .218 .774 1.344 1.566 2.119 ];
%   for n=1:length(ztide)
%               
%                  plot(xl,[ztide(n) ztide(n)],'k:','linewidth',2);
%                  plot(xl,[ztide(n) ztide(n)],'b:','linewidth',1); 
%   end
%         plot(xl,[ztide(2)-.25 ztide(2)-.25],'r--');
%         text(xl(end),ztide(1),' LAT','fontsize',14,'verticalalign','bottom');
%         text(xl(end),ztide(2)-0.25,' P','color','r','verticalalign','bottom');
%         text(xl(end),ztide(2),' MLLW','verticalalign','bottom');
%         text(xl(end),ztide(3),' MLW','verticalalign','bottom');
%         text(xl(end),ztide(4),' MSL','fontsize',14,'verticalalign','bottom');
%         text(xl(end),ztide(5),' MHW','verticalalign','bottom');
%         text(xl(end),ztide(6),' MHHW','verticalalign','bottom');
%         text(xl(end),ztide(7),' HAT','fontsize',14,'verticalalign','bottom');
%         
%         plot(SM(end).Xshoreline,ShoreElev,'k.','markersize',20)
%         
% %  legend([pl pe],datestr(SM(end).Datenum,'mm/dd/yyyy'),'Most Eroded',...
% %      'location','northwest')
%  
%  xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88');
%  
%  title({['Historical ' num2str(MopNumber) ' Profiles with Shoreline Location'];...
%      ['+/-' num2str(dx) 'm from most recent Survey']}) 
% 
% 
% legend([pl pe pc],datestr(SM(end).Datenum,'mm/dd/yyyy'),'Most Eroded',...
%      'past ATV cobble/bedrock flags','location','northwest')
%  
% makepng(['MLLW' num2str(MopNumber) 'location' num2str(dx) '.png'])
% 
% % % make cobble atv sighting history plot
% % figure;plot(x1d(col),dt,'m.');datetick('y');
% % set(gca,'xdir','reverse');grid on;
% % xlabel('Xshore Distance (m)');
% % ylabel('Survey Date');
% % set(gca,'fontsize',14);
% % title(['Historical ' num2str(MopNumber) ' Profiles with ATV Cobble/Bedrock Sightings']);
%      
% %makepng(['ATVCobbleBedrockFlags' num2str(MopNumber) '.png'])
% 
