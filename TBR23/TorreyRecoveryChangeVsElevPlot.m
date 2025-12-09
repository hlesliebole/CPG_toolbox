clearvars
%close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% loop throu TP mops 578-595
%close all
figure('position',[1          55        1412         742]);
hold on
m=0;
col=jet(10);
z1t=[];
z2t=[];
for MopNumber=580:589
    

load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));
%ndx=find([SA.Datenum] > datenum(2023,4,16));
%ndx=find([SA.Datenum] > datenum(2017,9,1) & [SA.Datenum] < datenum(2019,9,1));
% fall surveys
% mon=month(datetime([SA.Datenum],'convertfrom','datenum'));
% ndx=find(mon > 6 & mon < 10);
% stype=cell(numel(ndx),1);
% 
% all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};



%for n=ndx 
m=m+1;
n=jumbo(jstart);
n2=size(SA,2);n2=n+2;
     [x1d,z1]=GetNonGriddedProfile(MopNumber,n);
     [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
     
%      v(m)=sum(z1di(z1di > 2 & z1di < 3.5),'omitnan');
%      t(m)=SA(n).Datenum;
     %fprintf('%s %i\n')
%         if n == ndx(end)
%            
%             
%          p(m)=plot(x1d,z1di,'m*-',...
%         'DisplayName',[num2str(MopNumber) ' ' datestr(SA(n).Datenum)],...
%         'linewidth',4);hold on;
%     
%         elseif n == ndx(end-1)
           
%             %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
%          p(m)=plot(x1d,z1di,'m:',...
%         'DisplayName',[num2str(MopNumber) ' ' datestr(SA(n).Datenum)],...
%         'linewidth',4);hold on;
%         else
        %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
        p(m)=plot(z1-0.774,100*(z2-z1),'.','color',col(m,:),'markersize',10,...
         'DisplayName',[num2str(MopNumber)]);hold on;
        z1t=[z1t z1];z2t=[z2t z2];
%          p(m)=plot(x1d,movmean(z1di-z1di2,1,'omitnan'),'LineStyle','-',...
%              'color',col(m,:),...
%         'DisplayName',[num2str(MopNumber)],...
%         'linewidth',2);hold on;
        %end
%end
%   m=m+1;
%   [x1d,z1di]=GetNonGriddedProfile(MopNumber,2);
%   p(m)=plot(x1d,z1di,'k-','linewidth',4,'displayname','1998 El Nino');
% %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);

end

idx=find(z1t > -8 & z1t < -3);
bias=mean(z2t(idx)-z1t(idx),'omitnan');

% sort starting elevations
[z1s,i]=sort(z1t,'ascend');
% elev change of sorted starting elevations
zds=(z2t(i)-z1t(i));
%  smoothed change vs elevation
% 200 pts in sorted elevation is ~0.5m elevation change
pm=plot(z1s-0.774,movmean(100*zds,200,'omitnan'),'k-','linewidth',2,...
    'DisplayName','0.5m Moving Mean');

grid on;
box on;
set(gca,'linewidth',2);
set(gca,'fontsize',16);
xlabel(['Starting Elevation on ' datestr(SA(n).Datenum) ' (m, MSL)']);
ylabel('Elevation Change (cm)')
% set(gca,'ylim',[-.5 .5],'xdir','reverse','fontsize',16);
set(gca,'ylim',[-30 30],'ytick',-30:5:30,'xlim',[-11 -3]);
xl=get(gca,'xlim');plot(xl,[0 0],'k--','linewidth',2);
title([{'Mop 580-589 Average'},{['JetSki Survey Elev change since '...
    datestr(SA(n).Datenum) ' to ' datestr(SA(n2).Datenum) ]}],...
   'fontsize',18);
legend([fliplr(p) pm],'location','eastoutside','color',[.8 .8 .8])
set(gca,'color',[.8 .8 .8])
%text(-6.5,.22,['Bias (-8m < z < -3m) = ' num2str(bias,'%5.3fm')],'fontsize',18,'backgroundcolor','w')
xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
% plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
% plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
% plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',16,'linewidth',2);
% plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',16,'linewidth',2);  
%legend(p,'location','eastoutside','fontsize',14)
% legend(fliplr(p),'location','eastoutside','fontsize',14,'numcolumns',1);%ceil(size(SA,2)/60))
% set(gca,'position',[0.1143    0.5100    0.6817    0.4150])
% 
% axes
% PlotTideAndHs
% set(gca,'position',[0.1143    0.0550    0.6817    0.350])
% set(gca,'fontsize',14)
% box on;
% set(gca,'linewidth',2)
% yyaxis left
% xs=datetime(SA(end).Datenum,'convertfrom','datenum','timezone','utc');
% plot([xs xs],[-.5 2],'b-','linewidth',3);
% xs2=datetime(SA(end-1).Datenum,'convertfrom','datenum','timezone','utc');
% plot([xs2 xs2],[-.5 2],'b-','linewidth',3);
% fill([xs xs2 xs2 xs xs],[-.5 -.5 2 2 -.5],'y','FaceAlpha',.1)
% legend([pf pt],'location','eastoutside')
% set(gca,'position',[0.1143    0.0550    0.6817    0.350])
% makepng(['TorreyJetSkiProfileChange.png'])
%   
