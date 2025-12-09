addpath /Users/William/Desktop/Noah
MopNumber=582;
buoyid='100';
load CDIP100BP.mat
% fprintf('Building Buoy + Mop Prediction Data time series. Takes awhile...\n')
%     [Freq,Bw,TimeUTC,a0,a1,b1,a2,b2,Dtype]=getcdipBuoyAndModelnetcdf(buoyid);
% wavehs=4*sqrt(Bw'*a0); 
% wdate=TimeUTC;
% Get the deep water (set depth = 2000m here) group velocities (Cg)
%  of each of the buoy freqeuncies for event energy flux calculations
% depth=2000;
% [L,C,Cg]=LinearDispersion(Freq,2000);
% Cg=repmat(Cg,[1 size(a0,2)]); % convert group velcity vector in 2d array match a0events;
% Eflux=Bw'*(a0.*Cg); % total energy flux E(f)*Cg(f) in swell bands
% wavehs=Eflux;

%[wavetime,wavehs]=getwaves(MopNumber,datenum(2000,1,1),datenum(2023,1,1),'utc');
%wdate=datenum(wavetime);
% figure
% y=2008;
% idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
% subplot(2,1,1)
% plot(wavetime(idx),wavehs(idx))
% hold on
% y=2021;
% idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
% subplot(2,1,2)
% plot(wavetime(idx),wavehs(idx))

figure
legend
for y=2001:2022%[2008 2010 2014 2016 2017 2021]
    idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
    [n,edges]=histcounts(wavehs(idx),0:.1:6);
    plot(5.95:-0.1:0.05,cumsum(fliplr(n)),'linewidth',2,'DisplayName',num2str(y));...
        set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
    hold on;
end
for y=[2003 2007 2010 2016 2017 2021]
    idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
    [n,edges]=histcounts(wavehs(idx),0:.1:6);
    plot(5.95:-0.1:0.05,cumsum(fliplr(n)),'^-','linewidth',2,'DisplayName',num2str(y));...
        set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
    hold on;
end
  legend('location','northwest') 
  
figure
legend
for y=2001:2022%[2008 2010 2014 2016 2017 2021]
    idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
    [n,edges]=histcounts(Eflux(idx),0:.1:15);
    plot(14.95:-0.1:0.05,cumsum(fliplr(n)),'linewidth',2,'DisplayName',num2str(y));...
        set(gca,'xdir','reverse');set(gca,'ylim',[0 150]);grid on
    hold on;
end
for y=[2003 2007 2010 2016 2017 2021]
    idx=find(wdate > datenum(y-1,10,1) & wdate < datenum(y,10,1));
    [n,edges]=histcounts(Eflux(idx),0:.1:15);
    plot(14.95:-0.1:0.05,cumsum(fliplr(n)),'^-','linewidth',2,'DisplayName',num2str(y));...
        set(gca,'xdir','normal');set(gca,'ylim',[0 150]);grid on
    hold on;
end
  legend('location','northwest') 
% yr=year(wavetime);
% [n,edges]=histcounts(wavehs(yr == 2021),0:.1:6);
% 
% [n,edges]=histcounts(wavehs(yr == 2020),0:.1:6);
% plot(5.95:-0.1:0.05,cumsum(fliplr(n)));set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on
% hold on;
% [n,edges]=histcounts(wavehs(yr == 2003),0:.1:6);
% plot(5.95:-0.1:0.05,cumsum(fliplr(n)));set(gca,'xdir','reverse');set(gca,'ylim',[0 100]);grid on

