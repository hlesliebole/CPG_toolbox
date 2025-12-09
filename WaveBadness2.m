addpath /Users/William/Desktop/Noah
MopNumber=582;
buoyid='100';
%load CDIP100BP.mat
% fprintf('Building Buoy + Mop Prediction Data time series. Takes awhile...\n')
%     [Freq,Bw,TimeUTC,a0,a1,b1,a2,b2,Dtype]=getcdipBuoyAndModelnetcdf(buoyid);
% save TorreyOuterFilled.mat Freq Bw TimeUTC a0 a1 b1 a2 b2 Dtype    
load TorreyOuterFilled.mat    
wavehs=4*sqrt(Bw'*a0); 
wdate=TimeUTC;
% Get the deep water (set depth = 2000m here) group velocities (Cg)
%  of each of the buoy freqeuncies for event energy flux calculations
depth=2000;
[L,C,Cg]=LinearDispersion(Freq,2000);
Cg=repmat(Cg,[1 size(a0,2)]); % convert group velcity vector in 2d array match a0events;
Eflux=Bw'*(a0.*Cg); % total energy flux E(f)*Cg(f) in swell bands
close all

%figure;


minHs=2.0;
dh=0.5;
maxHs=max(wavehs);
maxHs=dh*(ceil(maxHs/dh));

figure('position',[ 317   223   927   554]);
legend
for y=2001:2023%[2008 2010 2014 2016 2017 2021]
    m=0;
    for hs=minHs:dh:maxHs
     m=m+1;
     d(m)=hs*2;
     idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1)...
         & wavehs > hs & wavehs <= hs+dh);
     Efsum(m)=sum(Eflux(idx),'omitnan');
    end
    
    %[n,edges]=histcounts(wavehs(idx),0:.1:6);
    if y == 2016
     plot(d,Efsum,'m.-','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' El Nino Year']);...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    elseif y == 2017
     plot(d,Efsum,'m.:','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' La Nina Year']);...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    elseif y == 2021
     plot(d,Efsum,'b.-','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' Extreme Event Year']);...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    elseif y == 2010
     plot(d,Efsum,'g.-','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' El Nino Year']);...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    elseif y == 2023
     plot(d,Efsum,'r.:','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' This Year']);...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
%     elseif y == 2002
%      plot(d,Efsum,'c:','linewidth',3,'markersize',15,'DisplayName',['\bf' num2str(y) ' El Nino Year']);...
%         set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    else 
     plot(d,Efsum,'.-','color',[.8 .8 .8],'linewidth',2,'markersize',15,'DisplayName',num2str(y));...
        set(gca,'xdir','normal');%set(gca,'ylim',[0 100]);grid on
    end
    hold on;
end
grid on;
set(gca,'fontsize',16);
title('San Diego County Local Deepwater Waves')
xlabel('\bfDeposition Depth (m)','fontsize',18);
ylabel([{'\bfWave Deposition Potential'},{'\fontsize{14}Net Hourly Wave Power (m^{3}/s) '}],'fontsize',18 )

% for y=[2003 2007 2010 2016 2017 2021]
%     idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
%     [n,edges]=histcounts(wavehs(idx),0:.1:6);
%     plot(5.95:-0.1:0.05,cumsum(fliplr(n)),'^-','linewidth',2,'DisplayName',num2str(y));...
%         set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
%     hold on;
% end
%   legend('location','northwest') 
%   
% figure
% legend
% for y=2001:2022%[2008 2010 2014 2016 2017 2021]
%     idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
%     [n,edges]=histcounts(Eflux(idx),0:.1:15);
%     plot(14.95:-0.1:0.05,cumsum(fliplr(n)),'linewidth',2,'DisplayName',num2str(y));...
%         set(gca,'xdir','reverse');set(gca,'ylim',[0 150]);grid on
%     hold on;
% end
% for y=[2003 2007 2010 2016 2017 2021]
%     idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
%     [n,edges]=histcounts(Eflux(idx),0:.1:15);
%     plot(14.95:-0.1:0.05,cumsum(fliplr(n)),'^-','linewidth',2,'DisplayName',num2str(y));...
%         set(gca,'xdir','normal');set(gca,'ylim',[0 150]);grid on
%     hold on;
% end
  legend('location','northeast') 
% yr=year(wavetime);
% [n,edges]=histcounts(wavehs(yr == 2021),0:.1:6);
% 
% [n,edges]=histcounts(wavehs(yr == 2020),0:.1:6);
% plot(5.95:-0.1:0.05,cumsum(fliplr(n)));set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
% hold on;
% [n,edges]=histcounts(wavehs(yr == 2003),0:.1:6);
% plot(5.95:-0.1:0.05,cumsum(fliplr(n)));set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
makepng('DepositionPotential.png')
