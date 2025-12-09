% Example code to make ENSO pattern annual mean MHW beach width plot
%  for a specified mop range compared to Sat estimate
%   currently only works between Mops 497-682 (only range I preprocessed)
%   Coast Sat data

%% uncomment the desired beach to plot

MopRange=[666 682];ttl='Cardiff SB';fn='CardiffvsVos.png';
% MopRange=[639 663];ttl='Solana Beach';fn='SolanavsVos.png';
% MopRange=[621 635];ttl='North Del Mar';fn='NorthDelMarvsVos.png';
% MopRange=[591 620];ttl='South Del Mar';fn='SouthDelMarvsVos.png';
% MopRange=[579 590];ttl='North Torrey';fn='NorthTorreyvsVos.png';
% MopRange=[540 565];ttl='South Torrey';fn='SouthTorryvsVos.png';
% MopRange=[497 514];ttl='La Jolla Shores';fn='LJShoresvsVos.png';


dname={[ttl ' : MOPs ' num2str(MopRange(1)) ' to '  num2str(MopRange(2))],...
    'Annual Mean MHW Beach Width'};

% call function to get annual widths for this range
[xmbw,mbw]=GetReachAnnualMeanBeachWidths(MopRange);

% make a plot
figure('position',[283   360   846   437]);

% xmbw are beach years, offset value by .3 years to center on oct-sep
%  beach year time period
h=plot(xmbw+.3,mbw,'k*','linewidth',2,'markersize',10,'DisplayName',...
    'SIO Survey-derived');hold on;
hold on;

% plot linear fit trend lines between EL Nino years using
%   function linfit
j=find(xmbw > 2015 & xmbw < 2021); % end before 2021 extreme event year
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
   col=get(h,'color'); 
   xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end
    j=find(xmbw > 1997 & xmbw < 2003);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end  
    j=find(xmbw > 2002 & xmbw < 2010);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end  
    j=find(xmbw > 2009 & xmbw < 2016);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));lf=plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2,...
        'DisplayName','SIO Survey Trend');
end  


xlabel('Beach Year')
ylabel('MHW Beach Width (m)');
%set(gca,'ylim',yl);
set(gca,'xlim',[1997 2023],'xtick',1997:2023,'fontsize',14,'linewidth',2);
grid on

title(dname,'fontsize',20)
%makepng('Test.png')

xsio=xmbw;
bwsio=mbw;
% call function to get annual widths for this range
[xmbw,mbw]=GetReachAnnualMeanBeachWidthsVos(MopRange);

h2=plot(xmbw+.3,mbw,'ro','linewidth',2,'markersize',10,'DisplayName',...
    'Satellite-derived');hold on;


% gray shade el nino years 
ey=1998;

% plot linear fit trend lines between EL Nino years using
%   function linfit
col='r';
j=find(xmbw > 2015 & xmbw < 2021); % end before 2021 extreme event year
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
   xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end
    j=find(xmbw > 1997 & xmbw < 2003);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end  
    j=find(xmbw > 2002 & xmbw < 2010);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end  
    j=find(xmbw > 2009 & xmbw < 2016);
if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
    xy=min(xmbw(j)):max(xmbw(j));lf2=plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2,...
        'DisplayName','Satellite Trend');
end  

yl=get(gca,'ylim');
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2003;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2010;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2016;
en=fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8],...
    'DisplayName','El Nino Year');


legend([h lf h2 lf2 en],'location','southoutside','numcolumns',3)

[val,ia,ib]=intersect(xmbw,xsio);
bias=mean(mbw(ia)-bwsio(ib),'omitnan');
text(2004,yl(1)-(yl(2)-yl(1))*.6,['Mean Satellite minus Survey Difference: ' num2str(bias,'%4.1fm')],'backgroundcolor','w','fontsize',16 )
%makepng('CardiffBeachWidthSIOvsVos.png')
makepng(fn)