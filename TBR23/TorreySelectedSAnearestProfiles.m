clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% plot Mop average profiles for slected dates to illustrate
%   Torrey Pines Recovery
Jdates=[datenum(2023,4,6) datenum(2023,8,15) datenum(2023,10,16) datenum(2023,12,12) datenum(2024,1,10) datenum(2024,1,26)];

% figure out how many surveys to process starting 6 Apr
MopNumber=584;
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski') | contains({SA.Source},'Trk' ));
[X1D,Z1D]=GetNearestPointsProfiles(SA,25,5);

nsurveys=numel(Jdates);
col=jet(nsurveys);
% col(1,:)=[0 0.7 0];
% col(2,:)=[0 0 0 ];
% col(3,:)=[0.5 0.5 0.9 ];
% col(4,:)=[1 0 0];
% col(5,:)=[.7 .7 .7];


figure('position',[ 42          85        1307         710])
for n=1:nsurveys
stn=find([SA(jumbo).Datenum] == Jdates(n));
for sn=stn
x1d=X1D;
z1=Z1D(jumbo(sn),:);

%[x1d,z1]=GetNonGriddedProfile(MopNumber,jumbo(sn));
% remove 20cm bias from Apr 2023 survey
if Jdates(n) == datenum(2023,4,6)
    z1=z1-.2;
end
p(n)=plot(x1d,movmean(z1-0.774,7,'omitnan'),'-','color',col(n,:),'linewidth',3,'DisplayName',datestr(Jdates(n)));
hold on
end
end
plot([0 620],[0 0],'k-')
% plot([220 620],[-4 -4],'k:')
% plot([220 220],[-12 -4],'k:')
set(gca,'xdir','reverse','xlim',[0 620],'fontsize',18,'color',[.8 .8 .8]);grid on;
legend(p,'location','north','fontsize',20)
title(['Torrey Pines 2023-24 : Mop ' num2str(MopNumber)],'fontsize',22)
%title({'Torrey Pines 2022-23','Severe Winter Erosion and Partial Summer Recovery'},'fontsize',22);
xlabel('Cross-Shore Distance (m)');ylabel('Elevation (m, MSL)');
set(gcf,'inverthardcopy','off');
makepng(['TorreyMop' num2str(MopNumber) 'Profiles.png'])
% 
% fill([220 620 620 220 220],[-12 -12 -4 -4 -12],'y','facealpha',.2)
% legend(p,'location','north','fontsize',20)
% 
% text(570,-3.5,'Weak Recovery Below -4m MSL','fontsize',22)
% makepng('Torrey2023ProfilesHighlight.png')



% 
% % loop throu TP mops 578-595
% %close all
% figure('position',[1          55        1412         742]);
% hold on
% m=0;
% 
% for mm=1:nsurveys+3%nsurveys+1
%     m=m+1;
% z1t=[];
% z2t=[];    
% for MopNumber=580:589%589
% 
% load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
% jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
% jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));
% 
% % get starting survey
% n=jumbo(jstart);
% [x1d,z1]=GetNonGriddedProfile(MopNumber,n);
% z1=z1-.2;% assume -20 cm bias in initial surveys
% 
% if mm == nsurveys+1
%     n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,4,15)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2023,2,10)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2017,4,11)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,7,11)));
%     [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
% elseif mm == nsurveys+2
%     n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,7,11)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2023,2,10)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2023,1,24)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2017,7,19)));
%     % n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,10,10)));
%     [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
% elseif mm == nsurveys+3
%     n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,10,10)));
%     %n2=jumbo(find([SA(jumbo).Datenum] == datenum(2017,10,4)));
%     [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
% else
% % get next sequential survey after that
% n2=jumbo(jstart+m);
% [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
% end     
% %         p(m)=plot(z1-0.774,100*(z2-z1),'.','color',col(m,:),'markersize',10,...
% %          'DisplayName',[num2str(MopNumber)]);hold on;
% 
% % add to totals
% z1t=[z1t z1];z2t=[z2t z2];
% 
% end
% 
% idx=find(z1t > -8 & z1t < -3);
% bias=mean(z2t(idx)-z1t(idx),'omitnan');
% 
% % sort starting elevations
% [z1s,i]=sort(z1t,'ascend');
% % elev change of sorted starting elevations
% zds=(z2t(i)-z1t(i));
% %numel(find(z1s > -8 & z1s < -7.5 & ~isnan(zds)))
% zmm=movmean(zds,120,'omitnan');
% ibias=find(z1s-0.774 > -10 & z1s-0.774 < -9 & ~isnan(zmm));
% bias=-mean(zmm(ibias));
% zmm=zmm+bias;
% %  smoothed change vs elevation
% % 120 pts in sorted elevation is ~0.5m elevation change
% if mm == nsurveys+1
%     pm(m)=plot(z1s-0.774,movmean(100*zds,120,'omitnan'),'m-','linewidth',2,...
%     'DisplayName',[datestr(SA(n2).Datenum)]);
% hold on;
% elseif mm == nsurveys+2
%     pm(m)=plot(z1s-0.774,movmean(100*zds,120,'omitnan'),'k-','linewidth',2,...
%     'DisplayName',[datestr(SA(n2).Datenum)]);
% hold on;
% elseif mm == nsurveys+3
%     pm(m)=plot(z1s-0.774,movmean(100*zds,120,'omitnan'),'k:','linewidth',2,...
%     'DisplayName',[datestr(SA(n2).Datenum)]);
% hold on;
% else
% pm(m)=plot(z1s-0.774,zmm*100,'-','linewidth',2,...
%     'color',col(m,:),'DisplayName',[datestr(SA(n2).Datenum) ' Bias: ' num2str(-bias*100,'%4.1fcm')]);
% hold on;
% end
% end
% 
% grid on;
% box on;
% set(gca,'linewidth',2);
% set(gca,'fontsize',16);
% xlabel(['Starting Elevation on ' datestr(SA(n).Datenum) ' (m, MSL)']);
% ylabel('Elevation Change (cm)')
% % set(gca,'ylim',[-.5 .5],'xdir','reverse','fontsize',16);
% set(gca,'ylim',[-140 90],'ytick',-140:10:90,'xlim',[-11 -3]);
% xl=get(gca,'xlim');plot(xl,[0 0],'k--','linewidth',2);
% title([{'Mop 580-589 Average'},{['JetSki Survey Elev difference relative '...
%     datestr(SA(n).Datenum)]}],...
%    'fontsize',18);
% legend([pm],'location','eastoutside','color',[.8 .8 .8])
% set(gca,'color',[.8 .8 .8])
% set(gca,'position',[0.1143    0.5100    0.6817    0.4150])
% %text(-6.5,.22,['Bias (-8m < z < -3m) = ' num2str(bias,'%5.3fm')],'fontsize',18,'backgroundcolor','w')
% % text(-10.8,50,['Vertical Datum for each Survey is Bias Corrected to have 0 Mean '...Elev Change betwen 
% % 'Elevation Change Between -10m and -9m Depth.'],'fontsize',16,'backgroundcolor','w');
% 
% axes('position',[0.1143    0.0550    0.6817    0.350]);
% PlotTideAndHs
% set(gca,'position',[0.1143    0.0550    0.6817    0.350])
% set(gca,'fontsize',14)
% box on;
% set(gca,'linewidth',2)
% yyaxis left
% xs=datetime(SA(jumbo(jstart)).Datenum,'convertfrom','datenum','timezone','utc');
% plot([xs xs],[-.5 2],'k-','linewidth',3)
% for ns=1:nsurveys
% n1=jumbo(jstart+ns-1);n2=jumbo(jstart+ns);
% xs=datetime(SA(n1).Datenum,'convertfrom','datenum','timezone','utc');
% xs2=datetime(SA(n2).Datenum,'convertfrom','datenum','timezone','utc');
% plot([xs2 xs2],[-.5 2],'-','linewidth',3,'color',col(ns,:));
% fill([xs xs2 xs2 xs xs],[-.7 -.7 2.2 2.2 -.7],col(ns,:),'FaceAlpha',.1,'edgecolor','none')
% end
% legend([pf pt],'location','eastoutside')
% set(gca,'position',[0.1143    0.0550    0.6817    0.350])
% set(gcf,'inverthardcopy','off');
% %makepng('TorreyRecoveryJetSkiChangeWith2022.png')
% 
