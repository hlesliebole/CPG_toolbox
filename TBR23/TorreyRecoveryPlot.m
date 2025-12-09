clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% loop throu TP mops 578-595
close all
figure('position',[1          55        1412         742]);
hold on
m=0;
col=jet(18);
for MopNumber=578:595
    

load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
ndx=find([SA.Datenum] > datenum(2023,4,18));
%ndx=find([SA.Datenum] > datenum(2017,9,1) & [SA.Datenum] < datenum(2019,9,1));
% fall surveys
% mon=month(datetime([SA.Datenum],'convertfrom','datenum'));
% ndx=find(mon > 6 & mon < 10);
stype=cell(numel(ndx),1);

all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};



for n=ndx 
m=m+1;     
     [x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
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
         p(m)=plot(x1d,z1di,'LineStyle','-',...
             'color',col(m,:),...
        'DisplayName',[num2str(MopNumber) ' ' datestr(SA(n).Datenum)],...
        'linewidth',2);hold on;
        %end
end
%   m=m+1;
%   [x1d,z1di]=GetNonGriddedProfile(MopNumber,2);
%   p(m)=plot(x1d,z1di,'k-','linewidth',4,'displayname','1998 El Nino');
% %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);

end

load M00587SM.mat;
jumbo=find(contains({SM.File},'umbo'));
m=m+1;
hold on;p(m)=plot(SM(jumbo(end-4)).X1D,SM(jumbo(end-4)).Z1Dtransect',...
    'm:','linewidth',3,...
    'DisplayName',[num2str(587) ' ' datestr(SM(jumbo(end-4)).Datenum)]);
m=m+1;
hold on;p(m)=plot(SM(jumbo(end-3)).X1D,SM(jumbo(end-3)).Z1Dtransect',...
    'k:','linewidth',3,...
    'DisplayName',[num2str(587) ' ' datestr(SM(jumbo(end-3)).Datenum)]);
grid on;
box on;
set(gca,'linewidth',2);
set(gca,'fontsize',16);xlabel('Distance From Mop Back Beach Point (m)');ylabel('Elevation (m, NAVD88)')
set(gca,'ylim',[-9 1.5],'xdir','reverse','xlim',[100 600],'fontsize',16);

title([{'Mops 578-595'},{['JetSki Survey on ' datestr(SA(ndx(1)).Datenum)]}],...
   'fontsize',18);
xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
% plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
% plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',16,'linewidth',2);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',16,'linewidth',2);  
%legend(p,'location','eastoutside','fontsize',14)
legend(fliplr(p),'location','eastoutside','fontsize',14,'numcolumns',1);%ceil(size(SA,2)/60))
set(gca,'position',[0.1143    0.5100    0.6817    0.4150])

axes
PlotTideAndHs
set(gca,'position',[0.1143    0.0550    0.6817    0.350])
set(gca,'fontsize',14)
box on;
set(gca,'linewidth',2)
yyaxis left
xs=datetime(SA(end).Datenum,'convertfrom','datenum','timezone','utc');
plot([xs xs],[-.5 2],'b-','linewidth',3);
legend([pf pt],'location','eastoutside')
set(gca,'position',[0.1143    0.0550    0.6817    0.350])
makepng(['TorreyJetSkiProfiles.png'])
  
