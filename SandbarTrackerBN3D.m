clear all
close all
mpath = '/volumes/group/MOPS/';
addpath(mpath);
addpath(fullfile(mpath, 'toolbox'));
addpath(fullfile(mpath, 'toolbox/profiles'));

MopNumber=582;%553;
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM')

% define global mean and global trend line
zmin=-10;
zmax=5;
SmoothWindow=11;
ijumbo=find(contains({SM.File},'umbo') | contains({SM.File},'jetski') ); % find jumbo survey
iwhl=find(contains({SM.Source},'wheel')); % find wheelie surveys
Zmean=mean(vertcat(SM(ijumbo).Z1Dtransect),'omitnan'); % global mean jumbo profile
X=SM(1).X1D; % xshore distance
Zs=Zmean( Zmean > zmin & Zmean < zmax ); % fixed section of global mean profile
Xs=SM(1).X1D( Zmean > zmin & Zmean < zmax ); % fixed section xshore distances
Zsdt=detrend(Zs,'omitnan'); % detrended global mean profile
Zsdtmm=movmean(Zsdt,SmoothWindow); % moving mean smoothing of detrended profile
GZsTrend=Zs-Zsdt; % global mean profile linear trend

bmax=1.5;bmin=-1.5;nc=16;
brange=bmax-bmin; % set max slope for coloring
cm =flipud(polarmap(nc));

GcsaX=[];
GcsaY=[];
VsaqX=[];
VsaX=[];
figure('position',[440   123   650   674]);hold on;
set(0, 'DefaultLineLineWidth', 2);
% loop through all the jumbo surveys 
for njumbo=1:numel(ijumbo)
    Z=SM(ijumbo(njumbo)).Z1Dtransect;
    y=SM(ijumbo(njumbo)).Datenum*ones(size(Z));
    ix=find(ismember(X,Xs));
    Zs=Z(ix);
    Zsdt=Zs-GZsTrend;
    %Zsdtmm=movmean(Zsaqdt,8);
    
    %m=Zsdt-min(Zsdt);
    m=Zsdt;
    %idx=find(~isnan(m) & GZsTrend > 0);
    idx=find(~isnan(m) & m > 0 & GZsTrend < 0);
    idxsa=find(~isnan(m) & GZsTrend > 0 & GZsTrend < 3);
    idxsaq=find(~isnan(m) & GZsTrend< 0);
    
    if ~isempty(idxsaq)
    Vsaq=sum(m(idxsaq));
    else
    Vsaq=NaN;
    end
    if ~isempty(idxsa)
    Vsa=sum(m(idxsa));
    else
    Vsa=NaN;
    end
    
    if ~isempty(idx)
    Xc=sum(m(idx).*Xs(idx))/sum(m(idx));
    else
    Xc=NaN;
    end
    [xmin,ixmin]=min(abs(Xs-Xc));
    Gcsa=GZsTrend(ixmin);
    
    idx=find(~isnan(m) & m > 0 );
    if ~isempty(idx)
    Xc=sum(m(idx).*Xs(idx))/sum(m(idx));
    else
    Xc=NaN;
    end
    [xmin,ixmin]=min(abs(Xs-Xc));
    Gcp=GZsTrend(ixmin);
    
%     idx=find(~isnan(m) & m < 0);
%     Xc=sum(m(idx).*Xs(idx))/sum(m(idx));
%     [xmin,ixmin]=min(abs(Xs-Xc));
%     Gcn=GZsTrend(ixmin);
    
    zscaled = 1+nc*(Zsdt-bmin)/brange;
    zscaled(zscaled < 1)=1;zscaled(zscaled > nc)=nc;
    idx=find(~isnan(zscaled)); % non NaN points 
    scp=scatter(GZsTrend(idx), y(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
    if min(Z) < -5
    p1=plot(Gcsa,y(1),'m.','markersize',15);%plot(Gcn,y(1),'mo');
    GcsaX=[GcsaX Gcsa];
    GcsaY=[GcsaY y(1)];
    VsaX=[VsaX Vsa];
    VsaqX=[VsaqX Vsaq];
    p2=plot(Gcp,y(1),'ko','markersize',5);
    end
    %fprintf('%6.1f %6.1f\n',min(Z),min(SM(ijumbo(njumbo)).Z1Dmean))
end

%plot(GcsaX,GcsaY,'m-')
yl=get(gca,'ylim');
datetick('y','mmm-YYYY');

ylabel('Survey Date');

%set(gca,'xdir','reverse','fontsize',12); 
%set(gca,'xlim',[min(Xs) max(Xs)]);
%x0=max(Xs(GZsTrend > 0));hold on;
x0=0;
plot([x0 x0],yl,'k-','linewidth',1);
set(gca,'xtick',[min(round(GZsTrend)):max(round(GZsTrend))],...
    'fontsize',12);
grid on %set(gca,'xgrid','on')
xlabel('Global Profile Linear Trend Line Elevation (m,NAVD88)')

cb=colorbar;
cb.Label.String='Detrended Elevation (m)';
cb.FontSize=12;
colormap(flipud(polarmap(nc)));
set(gca,'clim',[bmin bmax]);
title([{['Mop ' num2str(MopNumber) ' Detrended Profiles']},...
    {'(all relative to the global mean linear trend)'}],...
    'fontsize',16);
box on;


text(1,datenum(2022,10,1),'SUBAERIAL','fontsize',12,'fontweight','bold')
text(-5,datenum(2022,10,1),'SUBAQUEOUS','fontsize',12,'fontweight','bold')

pos=get(gca,'position');
set(gca,'position',[pos(1) pos(2)+.06 pos(3) pos(4)-.06])
pos=get(gca,'position');
ax2=axes('position',[pos(1) pos(2)-.07 pos(3) 0]);
set(ax2,'xlim',[min(Xs) max(Xs)],'xdir','reverse','fontsize',12);
xlabel('Xshore Distance (m)');
legend(ax2,[p2 p1],'Positive Elev Center of Mass Entire Profile',...
    'Positive Elev Center of Mass Subaqueous Profile','location','southoutside','numcolumns',2);

pngfile=['Mop' num2str(MopNumber) 'DetrendedJumbos.png'];
makepng(pngfile)

figure;
plot(GcsaY,GcsaX,'b.-','markersize',15);datetick('x','yyyy');grid on
hold on;
yyaxis right
plot(GcsaY,VsaX,'g.-','markersize',15);hold on;
plot(GcsaY,VsaqX,'r.-','markersize',15)
plot(GcsaY,VsaqX+VsaX,'k.-','markersize',15)
% figure('position',[138         187        1110         582]);
% set(0, 'DefaultLineLineWidth', 2);
% ax1=axes('position',[0.1300 0.1100 0.7750 0.3150]);
% %plot(X,Zmean,Xsaq,Zsaq,Xsaq,Zsaqdtmm,Xsaq,GZsaqTrend,'k:')
% 
% pg=plot(X,Zmean,'-','color','m');hold on;pt=plot(Xsaq,GZsaqTrend,'k:','linewidth',3)
% set(gca,'xdir','reverse','fontsize',12);grid on;
% x0=max(Xsaq(GZsaqTrend > 0));
% plot([x0 x0],yl,'k-','linewidth',1);hold on;
% legend([pg pt],['Global Mean Profile of (' num2str(numel(ijumbo)) ') Jumbos'],'Global Mean Profile Linear Trend',...
%     'location','northwest','fontsize',12)
% hold on;
% % [pks,ipks]=findpeaks(Zsaqdtmm);
% % for n=1:numel(pks)
% % plot([Xsaq(ipks(n)) Xsaq(ipks(n))],[Zsaq(ipks(n)) Zsaqdtmm(ipks(n))],...
% %     'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
% % end
% 
% xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
% title([' Global Mean Jumbo Profile'],'fontsize',16);
% set(gca,'xlim',[min(Xsaq(~isnan(Xsaq))) max(Xsaq(~isnan(Xsaq)))]);
% 
% 
% xl=get(gca,'xlim');
% yl=get(gca,'ylim');
% 
% ax2=axes('position',[0.1300 0.5100 0.7750 0.4550]);
% %plot(xl,[0 0],'k:');
% hold on;
% x0=max(Xsaq(GZsaqTrend > 0));
% plot([x0 x0],yl,'k-','linewidth',1);hold on;
% pgm=plot(Xsaq,Zsaqdtmm,'-','color','m','linewidth',3,...
%     'DisplayName','Global Mean Detrended Profile');%[0.8 0.8 0.8]);
% 
% for njumbo=1:numel(ijumbo)
%     Z=SM(ijumbo(njumbo)).Z1Dtransect;
%     ix=find(ismember(X,Xsaq));
%     Zsaq=Z(ix);
%     Zsaqdt=Zsaq-GZsaqTrend;
%     Zsaqdtmm=movmean(Zsaqdt,8);
%     mon=month(datetime(SM(ijumbo(njumbo)).Datenum,'convertfrom','datenum'));
%     yr=year(datetime(SM(ijumbo(njumbo)).Datenum,'convertfrom','datenum'));
%     if mon ==4 || mon == 5
%         if yr == 2011
%      ps=plot(Xsaq,Zsaqdtmm,'-','color','g','linewidth',3,...
%          'DisplayName',[datestr(SM(ijumbo(njumbo)).Datenum) ' Jumbo']);
%         end
%     else
%      ph=plot(Xsaq,Zsaqdtmm,'-','color',[.7 .7 .7],'linewidth',1,...
%          'DisplayName','Historic Jumbos 2004-2018');
%     end
%     if mon < 4
%         
%         if yr == 2006
%      pmx=plot(Xsaq,Zsaqdtmm,'-','color','r','linewidth',3,...
%          'DisplayName',[datestr(SM(ijumbo(njumbo)).Datenum) ' Badest Sandbar']);
%         axes(ax1)
%         plot(Xsaq,Zsaq,'-','color','r','linewidth',3,...
%          'DisplayName',[datestr(SM(ijumbo(njumbo)).Datenum) ' Badest Sandbar']);
%         axes(ax2)
%         end
%     end
% end
% pj=plot(Xsaq,Zsaqdtmm,'b-','linewidth',3,...
%     'DisplayName','4/20 Jetski & 4/19 Wheelie');
% 
% % find wheelie date closest to jetski date
% [dmin,imin]=min(abs([SM(iwhl).Datenum]-SM(ijumbo(njumbo)).Datenum)); 
% Z=SM(iwhl(imin)).Z1Dmean;
% ix=find(ismember(X,Xsaq));
% Zsaq=Z(ix);
% Zsaqdt=Zsaq-GZsaqTrend;
% Zsaqdtmm=movmean(Zsaqdt,8);
% plot(Xsaq,Zsaqdtmm,'b-','linewidth',3)
% 
% % most recent wheelie
% 
% %dmin,imin]=min(abs([SM(iwhl).Datenum]-SM(ijumbo(njumbo)).Datenum)); 
% Z=SM(iwhl(end)).Z1Dmean;
% ix=find(ismember(X,Xsaq));
% Zsaq=Z(ix);
% Zsaqdt=Zsaq-GZsaqTrend;
% Zsaqdtmm=movmean(Zsaqdt,8);
% plw=plot(Xsaq,Zsaqdtmm,'c-','linewidth',3,...
%     'DisplayName',[datestr(SM(iwhl(end)).Datenum) ' most recent wheelie']);
% 
% plot(xl,[0 0],'k:','linewidth',3);
% 
% text(100,-1.8,'SUBAERIAL','fontsize',14,'fontweight','bold')
% text(350,-1.8,'SUBAQUEOUS','fontsize',14,'fontweight','bold')
% legend([ph pgm pmx ps pj plw],'location','southwest');
% set(gca,'xdir','reverse','fontsize',12);grid on;box on
% set(gca,'xlim',xl);
% set(gca,'ylim',[-2 2]);
% ylabel('Detrended Elevation (m)');
% title(['Mop ' num2str(MopNumber) ' Detrended Profiles (all relative to the global mean linear trend)'],'fontsize',16);
% 
% pngfile=['Mop' num2str(MopNumber) 'GlobalMeanJumboBadestSandbar.png'];
% makepng(pngfile)
% 
% 
% % % process last jumbo using global trend line
% % %njumbo=find(ijumbo==24);
% % %njumbo=42; % Oct 2015
% % %njumbo=43; % Jan 2016
% % njumbo=numel(ijumbo); % latest jumbo
% % 
% % % find wheelie date closest to jetski date
% % [dmin,imin]=min(abs([SM(iwhl).Datenum]-SM(ijumbo(njumbo)).Datenum)); 
% % Z=SM(ijumbo(njumbo)).Z1Dtransect;
% % ix=find(ismember(X,Xsaq));
% % Zsaq=Z(ix);
% % Zsaqdt=Zsaq-GZsaqTrend;
% % Zsaqdtmm=movmean(Zsaqdt,8);
% % [pks,ipks,w,p]=findpeaks(Zsaqdtmm);
% % if ~isempty(pks) % delete small peaks relative the max peak prominence p
% %     ipks(p < 0.25*max(p))=[];
% %     pks(p < 0.25*max(p))=[];
% % end
% % plot(X,Z,'DisplayName',[datestr(SM(ijumbo(njumbo)).Datenum) ' Jetski Profile'])
% % plot(Xsaq,Zsaqdt,'DisplayName',[datestr(SM(ijumbo(njumbo)).Datenum) ' Detrended Profile'])
% % plot(X,SM(iwhl(imin)).Z1Dmean,'g-','DisplayName',[datestr(SM(iwhl(imin)).Datenum) ' Wheelie Profile'])
% 
% % figure('position',[138   122   613   647]);
% % set(0, 'DefaultLineLineWidth', 2);
% % plot(X,Z,Xsaq,Zsaqdt,Xsaq,GZsaqTrend,'k:')
% % set(gca,'xdir','reverse','fontsize',12);grid on;
% % legend(['Jumbo Profile ' datestr(SM(ijumbo(njumbo)).Datenum)],...
% %     'Detrended Jumbo Profile','Global Mean Subaqueous Profile Trend',...
% %     'location','southeast','fontsize',12)
% % hold on;
% % for n=1:numel(pks)
% % plot([Xsaq(ipks(n)) Xsaq(ipks(n))],[Zsaq(ipks(n)) Zsaqdtmm(ipks(n))],...
% %     'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
% % end
% % 
% % xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
% % title(['Mop ' num2str(MopNumber) ' Jumbo ' datestr(SM(ijumbo(njumbo)).Datenum)],'fontsize',16);
% % 
% % set(gca,'xlim',xl);
% % set(gca,'ylim',yl);
% % 
% % pngfile=['Mop' num2str(MopNumber) 'Jumbo' datestr(SM(ijumbo(njumbo)).Datenum,'YYYYmmDD') 'Sandbar.png'];
% % makepng(pngfile)
% % 
% % % 
% % % figure;plot(detrend(Zmean(Zmean < 0 & Zmean > -8),'omitnan'));hold on;
% % % plot(movmean(detrend(Zmean(Zmean < 0 & Zmean > -8)),11));
% % % plot(t,x,t,y,t,x-y,':k')
% % % legend('Input Data','Detrended Data','Trend','Location','northwest') 
% % % 
% % % figure;plot(detrend(SM(208).Z1Dmean,'omitnan'));hold on;plot(movmean(detrend(SM(208).Z1Dmean,'omitnan'),11))