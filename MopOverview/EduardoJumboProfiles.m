clearvars
close all

MopNumber=586;

% load the SA.mat file
SAmatfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
eval(['load ' SAmatfile]);

idx=find(contains({SA.File},'umbo') & strcmp({SA.Source},'Gps'));
SA=SA(idx);

NumSubTrans=1;
XgapTol=15;
YdistTol=40;
[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

% adjust to MSL
Z1Dtrans=Z1Dtrans-0.774;

col=jet(24);

figure('position',[30 61  1411 736]);

ax1=axes('position',[.05 .4 .4 .6]);
hold on;nn=0;
for syear=2001:2024
     nn=nn+1;
idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 4 |  month(Zdatetime) == 5));
%idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 12 |  month(Zdatetime) == 1));
ndx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 9 |  month(Zdatetime) == 10));
if numel(idx) > 0 & numel(ndx) > 0
idx=idx(ceil(numel(idx)/2));
ndx=ndx(ceil(numel(ndx)/2));
for n=idx   
    if year(Zdatetime(n)) == 2023
        plot(X1Dmop,movmean(Z1Dtrans(n,:),5,'omitnan'),'m-','linewidth',3,'DisplayName',datestr(SA(n).Datenum))
    else
        plot(X1Dmop,movmean(Z1Dtrans(n,:),5,'omitnan'),'-','color',col(nn,:),'linewidth',1,'DisplayName',datestr(SA(n).Datenum))
    end
end
end
end
set(gca,'xdir','reverse','fontsize',16);grid on;box on;
legend('location','northoutside','numcolumns',4);
ylabel('Elevation (m, MSL)');
xlabel('Cross-shore Distance (m)');
title(['MOP ' num2str(MopNumber) '  Spring Surveys'])


% -------------------------------------------------------

ax2=axes('position',[.55 .4 .4 .6]);
hold on;nn=0;
for syear=2001:2024
     nn=nn+1;
idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 9 |  month(Zdatetime) == 10));
ndx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 4 |  month(Zdatetime) == 5));
%ndx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 12 |  month(Zdatetime) == 1));
if numel(idx) > 0 & numel(ndx) > 0
idx=idx(ceil(numel(idx)/2));
ndx=ndx(ceil(numel(ndx)/2));
for n=idx   
    if year(Zdatetime(n)) == 2023
        plot(X1Dmop,movmean(Z1Dtrans(n,:),5,'omitnan'),'m-','linewidth',3,'DisplayName',datestr(SA(n).Datenum))
    else
        plot(X1Dmop,movmean(Z1Dtrans(n,:),5,'omitnan'),'-','color',col(nn,:),'linewidth',1,'DisplayName',datestr(SA(n).Datenum))
    end
end
end
end
set(gca,'xdir','reverse','fontsize',16);grid on;box on;
legend('location','northoutside','numcolumns',4)
title(['MOP ' num2str(MopNumber) '  Fall Surveys'])
ylabel('Elevation (m, MSL)')
xlabel('Cross-shore Distance (m)');

% -------------------------------------------------------

ax3=axes('position',[.05 .059 .4 .24]);
hold on;nn=0;
for syear=2001:2024
     nn=nn+1;
ndx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 9 |  month(Zdatetime) == 10));
idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 4 |  month(Zdatetime) == 5));
%idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 12 |  month(Zdatetime) == 1));
if numel(idx) > 0 & numel(ndx) > 0
idx=idx(ceil(numel(idx)/2));
ndx=ndx(ceil(numel(ndx)/2));
for n=1:numel(idx)   
    if year(Zdatetime(idx(n))) == 2023
        plot(X1Dmop,movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'m-','linewidth',3,'DisplayName',datestr(SA(n).Datenum))
    else
        plot(X1Dmop,movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'-','color',col(nn,:),'linewidth',1,'DisplayName',datestr(SA(n).Datenum))
    end
end
end
end
set(gca,'xdir','reverse','fontsize',16);grid on;box on;
%legend('location','northoutside','numcolumns',4)
title('Spring to Fall Elevation Change vs Cross-shore Distance')
ylabel('Change (m)')
xlabel('Cross-shore Distance (m)');

% -------------------------------------------------------

ax4=axes('position',[.55 .059 .4 .24]);
hold on;nn=0;
for syear=2001:2024
     nn=nn+1;
ndx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 9 |  month(Zdatetime) == 10));
idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 4 |  month(Zdatetime) == 5));
%idx=find(year(Zdatetime) == syear & ( month(Zdatetime) == 12 |  month(Zdatetime) == 1));
if numel(idx) > 0 & numel(ndx) > 0
idx=idx(ceil(numel(idx)/2));
ndx=ndx(ceil(numel(ndx)/2));
for n=1:numel(idx)  
    zmean=mean(vertcat(Z1Dtrans(idx(n),:),Z1Dtrans(ndx(n),:)));
    if year(Zdatetime(idx(n))) == 2023
        %plot(movmean(Z1Dtrans(idx(n),:),5,'omitnan'),movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'m-','linewidth',3,'DisplayName',datestr(SA(n).Datenum))
        plot(movmean(zmean,5,'omitnan'),movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'m-','linewidth',3,'DisplayName',datestr(SA(n).Datenum))
    else
        %plot(movmean(Z1Dtrans(idx(n),:),5,'omitnan'),movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'-','color',col(nn,:),'linewidth',1,'DisplayName',datestr(SA(n).Datenum))
        plot(movmean(zmean,5,'omitnan'),movmean(Z1Dtrans(ndx(n),:)-Z1Dtrans(idx(n),:),5,'omitnan'),'-','color',col(nn,:),'linewidth',1,'DisplayName',datestr(SA(n).Datenum))
    end
end
end
end
set(gca,'fontsize',16);grid on;box on;
%legend('location','northoutside','numcolumns',4)
title('Spring to Fall Elevation Change vs Mean Spring-Summer Depth')
ylabel('Change (m)')
xlabel('Mean Depth (m)');

exportgraphics(gcf,'TorreySpringFallProfileChange.png','Resolution',150)
