clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% figure out how many surveys to process starting 6 Apr
MopNumber=580;
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));
nsurveys=numel(jumbo)-jstart;
col=jet(nsurveys);

% loop throu TP mops 578-595
%close all
figure('position',[1          55        1412         742]);
hold on
m=0;

for mm=1:nsurveys%nsurveys+1
    m=m+1;
z1t=[];
z2t=[];    
for MopNumber=580:589%589
    
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));

% get starting survey
n=jumbo(jstart);
[x1d,z1]=GetNonGriddedProfile(MopNumber,n);

if mm > nsurveys
    n2=jumbo(find([SA(jumbo).Datenum] == datenum(2022,7,11)));
    [x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
else
% get next sequential survey after that
n2=jumbo(jstart+m);
[x1d,z2]=GetNonGriddedProfile(MopNumber,n2);
end     
%         p(m)=plot(z1-0.774,100*(z2-z1),'.','color',col(m,:),'markersize',10,...
%          'DisplayName',[num2str(MopNumber)]);hold on;

% add to totals
z1t=[z1t z1];z2t=[z2t z2];

end

idx=find(z1t > -8 & z1t < -3);
bias=mean(z2t(idx)-z1t(idx),'omitnan');

% sort starting elevations
[z1s,i]=sort(z1t,'ascend');
% elev change of sorted starting elevations
zds=(z2t(i)-z1t(i));
%numel(find(z1s > -8 & z1s < -7.5 & ~isnan(zds)))
zmm=movmean(zds,120,'omitnan');
ibias=find(z1s-0.774 > -10 & z1s-0.774 < -9 & ~isnan(zmm));
bias=-mean(zmm(ibias));
% Version with no bias correction
%zmm=zmm+bias;

%  smoothed change vs elevation
% 120 pts in sorted elevation is ~0.5m elevation change
if mm > nsurveys
    pm(m)=plot(z1s-0.774,movmean(100*zds,120,'omitnan'),'k-','linewidth',2,...
    'DisplayName',[datestr(SA(n2).Datenum) ' Bias:' num2str(-bias*100,'%4.1fcm')]);
hold on;
else
% pm(m)=plot(z1s-0.774,zmm*100,'-','linewidth',2,...
%     'color',col(m,:),'DisplayName',[datestr(SA(n2).Datenum) ' Bias: ' num2str(-bias*100,'%4.1fcm')]);
pm(m)=plot(z1s-0.774,zmm*100,'-','linewidth',2,...
    'color',col(m,:),'DisplayName',[datestr(SA(n2).Datenum) ]);
hold on;
end
end

grid on;
box on;
set(gca,'linewidth',2);
set(gca,'fontsize',16);
xlabel(['Starting Elevation on ' datestr(SA(n).Datenum) ' (m, MSL)']);
ylabel('Elevation Change (cm)')
% set(gca,'ylim',[-.5 .5],'xdir','reverse','fontsize',16);
set(gca,'ylim',[-45 20],'ytick',-40:5:15,'xlim',[-11 -3]);
xl=get(gca,'xlim');plot(xl,[0 0],'k--','linewidth',2);
title([{'Mop 580-589 Average'},{['JetSki Survey Elev change since '...
    datestr(SA(n).Datenum)]}],...
   'fontsize',18);
legend([pm],'location','eastoutside','color',[.8 .8 .8])
set(gca,'color',[.8 .8 .8])
set(gca,'position',[0.1143    0.5100    0.6817    0.4150])
%text(-6.5,.22,['Bias (-8m < z < -3m) = ' num2str(bias,'%5.3fm')],'fontsize',18,'backgroundcolor','w')
% text(-10.4,17,['Vertical Datum for each Survey is Bias Corrected to have 0 Mean '...Elev Change betwen 
% 'Elevation Change Between -10m and -9m Depth.'],'fontsize',16,'backgroundcolor','w');

axes('position',[0.1143    0.0550    0.6817    0.350]);
PlotTideAndHs
set(gca,'position',[0.1143    0.0550    0.6817    0.350])
set(gca,'fontsize',14)
box on;
set(gca,'linewidth',2)
yyaxis left
xs=datetime(SA(jumbo(jstart)).Datenum,'convertfrom','datenum','timezone','utc');
plot([xs xs],[-.5 2],'k-','linewidth',3)
for ns=1:nsurveys
n1=jumbo(jstart+ns-1);n2=jumbo(jstart+ns);
xs=datetime(SA(n1).Datenum,'convertfrom','datenum','timezone','utc');
xs2=datetime(SA(n2).Datenum,'convertfrom','datenum','timezone','utc');
plot([xs2 xs2],[-.5 2],'-','linewidth',3,'color',col(ns,:));
fill([xs xs2 xs2 xs xs],[-.5 -.5 2 2 -.5],col(ns,:),'FaceAlpha',.1,'edgecolor','none')
end
legend([pf pt],'location','eastoutside')
set(gca,'position',[0.1143    0.0550    0.6817    0.350])
set(gcf,'inverthardcopy','off');
makepng('TorreyRecoveryJetSkiChangeNoBiasCorrection.png')

