figure('position',[ 44         157        1261         613]);
dstart=datenum(2021,10,20);
dend=datenum(2022,1,25);
[tidetime,tidehgt]=getztide2(dstart,dend,'lst','msl','predictions');
xtide=datenum(tidetime);
plot(datetime(xtide,'convertfrom','datenum','TimeZone','America/Los_Angeles'),tidehgt,...
    '-','color',[.8 .8 .8],'linewidth',2);
ylabel('Tide (m, MSL)','color','k');set(gca,'ycolor','k');
set(gca,'fontsize',12)
hold on;
plot(datetime([dstart dend],'convertfrom','datenum','TimeZone','America/Los_Angeles'),[0 0],'k-');
ptide=tidehgt;ptide(tidehgt < 0)=NaN;
p1=plot(datetime(xtide,'convertfrom','datenum','TimeZone','America/Los_Angeles'),movmean(ptide,25,'omitnan'),...
    '-','color','g','linewidth',2);
ntide=tidehgt;ntide(tidehgt > 0)=NaN;
p2=plot(datetime(xtide,'convertfrom','datenum','TimeZone','America/Los_Angeles'),movmean(ntide,25,'omitnan'),...
    '-','color','r','linewidth',2);
legend([p1 p2],'25hr Moving Mean of Hourly Tide Elevations ABOVE MSL only',...
    '25hr Moving Mean of Hourly Tide Elevations BELOW MSL only','location','northwest')
xl=get(gca,'xlim');set(gca,'Xtick',xl(1):xl(2))
datetick('x','mm/dd','keepticks')
set(gca,'position',[.05 .1 .9 .8])
set(gca,'xlim',[datetime(min(xtide),'convertfrom','datenum','TimeZone','America/Los_Angeles')...
datetime(max(xtide),'convertfrom','datenum','TimeZone','America/Los_Angeles')]);
