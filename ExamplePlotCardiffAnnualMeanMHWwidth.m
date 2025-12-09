% Example code to make ENSO pattern annual mean MHW beach width plot
%  for a specified mop range

MopRange=[666 682]; % Cardiff SB Mop range
dname=['Cardiff Mops: ' num2str(MopRange(1)) ' to '  num2str(MopRange(2))];

% call function to get annual widths for this range
[xmbw,mbw]=GetReachAnnualMeanBeachWidths(MopRange);

% make a plot
figure('position',[283   426   846   291]);

% xmbw are beach years, offset value by .3 years to center on oct-sep
%  beach year time period
h=plot(xmbw+.3,mbw,'k*','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
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
    xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
end  

yl=get(gca,'ylim');
% gray shade el nino years 
ey=1998;

fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2003;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2010;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2016;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);

xlabel('Beach Year')
ylabel('Beach Width (m)');
set(gca,'ylim',yl);
set(gca,'xlim',[1997 2023],'xtick',1997:2023,'fontsize',14,'linewidth',2);
grid on
legend(h,'location','northwest')
title('Annual Mean MHW Beach Width Evolution','fontsize',22)
%makepng('Test.png')
