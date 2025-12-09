
% For an input global mean Mop profile elevation and Mop reach,
%  calculates the individual survey elevation anomalies at the xshore 
%  location(s) of the input elevation on each Mop's global mean profile.

%---------
% Settings
%---------
clear all



MopStart=520;  % mop range
MopEnd=596;
% MopStart=552;  % mop range
% MopEnd=552;
% running mean window (odd number of months) for mean seasonal MOP profiles
SeasonalAvgWindow=3;  

%----------------------

load MopTableUTM.mat

MopNumbers=MopStart:MopEnd;
NumMopPoints=length(MopNumbers);

GS.MopNumbers=MopNumbers;

% loop through Mop SM mat files

mm=0;
for m=MopStart:MopEnd  % Mop number loop
    mm=mm+1;
   
    % load mop SM mat file
    
    load([ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SM');
    load([ 'M'  num2str( m , '%5.5i' )  'GM.mat' ],'GM');
    
    fprintf('Loaded: %s %i\n',[ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],size(SM,2));
   
    dc=0.5; % contour spacing
    cnum=0;
 for ElvMsl=-9:dc:2   % elevation conotur loop
    cnum=cnum+1;
    %clear Anom
    
    %ElvMsl=2.5; % mean Mop profile elevation MSL for the anomaly estimates
     Elv=ElvMsl+0.774; % convert to navd88
    
    % get global mean and seaosnal mean profile Elv locations for this mop
       clear gx gsx
       
    % global Elv location(s) 
        xl=[GM.X1D(1) GM.X1D(end)];  
        %gx.x=intersections(xl,[Elv Elv],GM.X1D,GM.Z1Dmean);
        idx=find(GM.Z1Dmean >= Elv-0.5*dc & GM.Z1Dmean < Elv+0.5*dc);
        % round xshore locatiosn to nearest meter
        %gx.x=round(gx.x);
        if isempty(idx)
          gx.x=NaN;
        else
          gx.x=round(GM.X1D(idx));
        end
     
     % seasonal mean Elv location(s)
          
          for mon=1:12
              % make mean seasonal profile from monthly mean
              % profiles using SeasonalAvgWindow setting
              dm=(SeasonalAvgWindow-1)/2;
              mons=mon-dm:mon+dm;
              mons(mons < 1)=mons(mons < 1)+12;
              mons(mons > 12)=mons(mons > 12)-12;
              Zsm=nanmean(reshape([GM.MM(mons).Z1Dmean],[length(GM.MM(1).Z1Dmean),SeasonalAvgWindow]),2);
              xl=[GM.MM(mon).X1D(1) GM.MM(mon).X1D(end)];
              %dist=intersections(xl,[Elv Elv],GM.MM(mon).X1D,Zsm);
              %gsx.MM(mon).x=intersections(xl,[Elv Elv],GM.MM(mon).X1D,Zsm);
              idx=find(Zsm >= Elv-0.5*dc & Zsm < Elv+0.5*dc);
              %gsx.MM(mon).x=round(gsx.MM(mon).x);
              
%               if isempty(gsx.MM(mon).x)
%                   gsx.MM(mon).x=NaN;
%               end
              if isempty(idx)
                  gsx.MM(mon).x=NaN;
              else
                  gsx.MM(mon).x=GM.MM(mon).X1D(idx);
              end
          end
        
        
    % step through survey dates
    for n=1:size(SM,2)
        sdate=SM(n).Datenum;
        mon=month(datetime(datestr(sdate))); % survey month
        stype=SM(n).Source;
        
        % get anomaly relative to global mean profile
        idx=find(ismember(round(SM(n).X1D),gx.x));
        ganom=nanmean(SM(n).Z1Dmean(idx)-Elv)*length(idx);
        % get anomaly relative to global mean seasonal profile
        idx=find(ismember(round(SM(n).X1D),gsx.MM(mon).x));
        gsanom=nanmean(SM(n).Z1Dmean(idx)-Elv)*length(idx);
   
    %     % find index of exist Shoreline struct array entry for this date
    %  and survey source if it exists
        if m == MopStart && n == 1
          idx=[]; 
        else
         idx=find([Anom.Datenum] == sdate & strcmpi({Anom.Source}, stype)==1);
        end
    
    % if new date, add new entry to the Shoreline struct array and
    %  populate shoreline vectors with with NaNs for all Mops
      
       if isempty(idx)
           if m == MopStart && n == 1
           %if n == 1
               nn=1;
           else
               nn=size(Anom,2)+1;
           end
        Anom(nn).Datenum=sdate;
        Anom(nn).Source=stype;
%         Anom(nn).MopNumbers=MopStart:MopEnd;
%         Anom(nn).Global=nan(1,NumMopPoints);
%         Anom(nn).Seasonal=nan(1,NumMopPoints);
        
        Anom(nn).Global=ganom;
        Anom(nn).Seasonal=gsanom;
        Anom3D(nn,cnum,mm)=Anom(nn).Seasonal;
       else
        Anom(idx).Global=ganom;
        Anom(idx).Seasonal=gsanom;
        Anom3D(idx,cnum,mm)=Anom(idx).Seasonal;
       end
    end


% % add to Anom 3D (date,mop,contour,mop) array
% for dnum=1:size(Anom,2)
% Anom3D(dnum,cnum,mm)=Anom(dnum).Seasonal;
% end

end

end

d=[Anom.Datenum]; % survey dates in array
[sd,id]=sort(d);% sorted survey dates
jumbo=find(nansum(abs(squeeze(Anom3D(:,5,:))),2) > 0); % jumbo dates


%for ids=1:length(id)
nf=0;
clear M
figure('position',[233   301   784   477],'menu','none');
for ids=1:length(id)
%for ids=1:1
if(nansum(abs(Anom3D(id(ids),5,:))) > 0)
%close all
%figure('position',[233   301   784   377]);
nf=nf+1;
clf
ax1=axes('position',[.05 .87 .9 .05]);
set(ax1,'xlim',[min(d(jumbo)) max(d(jumbo))],'ylim',[-1 1]);
datetick;hold on;set(ax1,'ytick',[])
plot(d(jumbo),0*d(jumbo),'k.','markersize',10);hold on;
plot(sd(ids),0,'r.','markersize',20)
tx=text(sd(ids),1.8,datestr(sd(ids)),'horizontalalign',...
    'center','fontweight','bold','color','r');

ax2=axes('position',[.1 .3 .8 .5]);
colormap(flipud(polarmap));
a1=squeeze(Anom3D(id(ids),:,:));
imAlpha=ones(size(a1));
imAlpha(isnan(a1))=0;
%imAlpha(a1 < 0.1 )=0;
imagesc(520:596,-9:0.5:2,a1,'AlphaData',imAlpha);
hold on;plot([520 596],[0 0],'k--','linewidth',2);
c=colorbar;c.Label.String='Contour Xshore Volume Anomaly (m^3/m-shoreline)';
set(gca,'clim',[-25 25]);set(gca,'xticklabels',[])
set(gca,'Color',[.7 .7 .7]);grid on;
ylabel('Seasonal Mean Contour Elevation (m,MSL)');
title(datestr(sd(ids)))

ax3=axes('position',[.1 .1 .745 .18]);

set(gca,'xlim',[519.5 596.5],'ylim',[-310 310]);grid on;
hold on;plot([520 596],[0 0],'k--','linewidth',2);
ta=nansum(a1);pidx=find(ta > 0);nidx=find(ta < 0);
xt=520:596;
if ~isempty(pidx);plot(xt(pidx),ta(pidx),'b.','markersize',15);end
if ~isempty(nidx);plot(xt(nidx),ta(nidx),'r.','markersize',15);end
xlabel('Mop Number');
ylabel({'Total Xshore Vol','Anom (m^3/m-shore)'})
box on;
% estimate volume segments

if ~isempty(pidx)
tv=ta(pidx(1));
ip=pidx(1);
ps=pidx(1);
pl=pidx(1);
for i=2:length(pidx)
    if pidx(i)-ip == 1
        tv=nansum([tv ta(pidx(i))]);
        pl=pidx(i);
        ip=pidx(i);
    else
        if tv > 0 && (pl-ps) > 2
        if tv < 1
        text(xt(round((ps+pl)/2)),-200,...
         ['+' num2str(tv*100/1000,'%5.1f') 'K'],'fontweight','bold',...
         'color','b','horizontalalign','center');
        else
        text(xt(round((ps+pl)/2)),-200,...
        ['+' num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
        'color','b','horizontalalign','center'); 
        end
        end      
            tv=ta(pidx(i));
            ip=pidx(i);
            ps=pidx(i);
            pl=pidx(i);    
    end
end

if tv > 0 && (pl-ps) > 2
  if tv < 1
  text(xt(round((ps+pl)/2)),-200,...
   ['+' num2str(tv*100/1000,'%5.1f') 'K'],'fontweight','bold',...
   'color','b','horizontalalign','center');
  else
   text(xt(round((ps+pl)/2)),-200,...
   ['+' num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
   'color','b','horizontalalign','center');   
  end      
end
end
            
if ~isempty(nidx)
tv=ta(nidx(1));
ip=nidx(1);
ps=nidx(1);
pl=nidx(1);
for i=2:length(nidx)
    if nidx(i)-ip == 1
        tv=nansum([tv ta(nidx(i))]);
        pl=nidx(i);
        ip=nidx(i);
    else
        if tv < 0 && (pl-ps) > 2
        if tv > -1
        text(xt(round((ps+pl)/2)),200,...
         [num2str(tv*100/1000,'%5.1f') 'K'],'fontweight','bold',...
         'color','r','horizontalalign','center');
        else
        text(xt(round((ps+pl)/2)),200,...
        [num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
        'color','r','horizontalalign','center'); 
        end
        end 
            tv=ta(nidx(i));
            ip=nidx(i);
            ps=nidx(i);
            pl=nidx(i);    
    end
end
if tv < 0 && (pl-ps) > 2
  if tv > -1
  text(xt(round((ps+pl)/2)),200,...
   [num2str(tv*100/1000,'%5.1f') 'K'],'fontweight','bold',...
   'color','r','horizontalalign','center');
  else
   text(xt(round((ps+pl)/2)),200,...
   [num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
   'color','r','horizontalalign','center');   
  end      
end
end   

tv=nansum(ta);
text(600,200,[{'Net Volume','   Anomaly'}],'fontsize',12);
if tv < 0
  text(603,0,...
   [num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
   'color','r','verticalalign','middle',...
   'horizontalalign','center','fontsize',16);
  else
   text(603,0,...
   ['+' num2str(round(tv*100/1000)) 'K'],'fontweight','bold',...
   'color','b','verticalalign','middle',...
   'horizontalalign','center','fontsize',16);   
end  



M(nf)=getframe(gcf);
%pause
end
end
v=VideoWriter('TorreyAnomalySeasonal.mp4','MPEG-4');
v.FrameRate=1;
open(v)
writeVideo(v,M)
close(v)
% %[x y z]=ind2sub(size(Anom3D),find(~isnan(Anom3D)));
% %a=Anom3D(~isnan(Anom3D(:)));
% [x z y]=ind2sub(size(Anom3D),find(Anom3D > 0));
% a=Anom3D(Anom3D(:) > 0);
% % 
% % figure('position',[415   244   469   404],'name','Mop Mean Contour Anomalies');
% % x=repmat([Anom.Datenum],NumMopPoints,1);x=x(:)'; % make survey date for every points
% % y=[Anom.MopNumbers];
% % z=[Anom.Seasonal];
% % amax=ceil(max(abs(z)));
% amax=2;
% zscaled = 1+64*(a+amax)/(2*amax);
% zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
% idx=find(~isnan(zscaled)); % non NaN points
%                                      
% cm =flipud(polarmap(64));
% 
% scp=scatter3(x(idx), y(idx), z(idx), 10, cm(ceil(zscaled(idx)),:), 'filled');
% scp.MarkerFaceAlpha = .3;
% scp.MarkerEdgeAlpha = .3;
% %view(-10,80)
% colormap(flipud(polarmap(64)))
% cb=colorbar;cb.Label.String='Elevation Anomaly (m)';
% set(gca,'clim',[-amax amax]);
% set(gca,'xlim',[min(x) max(x)]);
% set(gca,'zlim',[min(z) max(z)]);
% set(gca,'color',[.7 .7 .7]);
% %datetick;
% xlabel('Survey Date');ylabel('Mop Number');zlabel('Anomaly (m)');
% title(['Elevation Anomalies at the ' num2str(ElvMsl) 'm MSL Seasonal Mean Contour']); 
% % view(2)
% % set(gcf,'inverthardcopy','off');
