% Example profile analysis focused on the the state of a Mop profile at
% time A relative to its monthly long-tmer mean and min,max,std

% close all
% clearvars

NumSubTrans=51;
XgapTol=5;
YdistTol=5;
XbackgapTol=15;

SBvol=NaN(31,numel(datenum(2024,1,1):datenum(2024,9,1)));
mn=0;
for MopNum=635:665
    MopNum
    mn=mn+1;

load(['M' num2str(MopNum,'%5.5i') 'SA.mat']);

% reduce to dates since Jan 1 2024
idx=find([SA.Datenum] >= datenum(2023,1,1) &...
    [SA.Datenum] <= datenum(2023,8,2) &...
    (strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR')) );
SA=SA(idx);

% toss 4 and 16 Apr at inlet Mop 637
if MopNum == 637
idx=find([SA.Datenum]==datenum(2023,4,17));SA(idx)=[];
idx=find([SA.Datenum]==datenum(2023,7,6));SA(idx)=[];
idx=find([SA.Datenum]==datenum(2023,8,1));SA(idx)=[];
end
% idx=find([SA.Datenum]==datenum(2024,2,12));SA(idx)=[];

% get profiles
[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

% [TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
%     TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%

Znavd88=0.744%1.344;
%[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,TZ1Dmon); 


[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,Z1Dmedian);
% figure;plot(Zdatetime,BeachVolMin,'+-',Zdatetime,BeachVolMax,'o-')
% title(' MHW Beach Vol');grid on;
for i=1:numel(Zdatetime)
    if ~isnan(BeachVolMax(i))
     iday=round(datenum(Zdatetime(i)))-datenum(2022,12,31);
     SBvol(mn,iday)=BeachVolMax(i);
    end
end

end

%%
SBvolMops=635:665;
SBvolDates=datenum(2023,1,1):datenum(2023,9,1);
figure
imagesc(SBvol','alphadata',~isnan(SBvol'));set(gca,'ydir','normal')
colormap(jet);
set(gca,'clim',[0 50]);

% vol matrix dates with data
igood=find(sum(SBvol,'omitnan') > 0);
SBdvol=SBvol*NaN;
SBdvol(:,igood(2:end))=diff(SBvol(:,igood),1,2);
%%
figure('position',[95         199        1090         598]);

col=flipud(jet(12));

ax1=axes('position',[0.075 0.05 .9 .4]);
for n=14:-1:2
    stn=['D' num2str(634+n,'%4.4i')];
        [wavetime,Sxy]=GetMopSxy(stn);
   
            
    if n == 2
            plot(wavetime,Sxy,...
        '-','color','m','linewidth',1,'displayname',['Mop ' num2str(634+n)])
        hold on;
    else
            plot(wavetime,Sxy,...
        '-','color',col(n-2,:),'linewidth',1,'displayname',['Mop ' num2str(634+n)])
        hold on;
    end
end
grid on;
set(gca,'xlim',[datetime(2023,1,1) datetime(2023,8,1)])
xtickformat('MMM-dd');set(gca,'fontsize',16,'color',[.8 .8 .8])
ylabel('Tot Sxy (m^2) + = North')
lg=legend('location','eastoutside');lg.Color=([ .8 .8 .8]);

ax2=axes('position',[0.075 0.55 .9 .4]);

stot=zeros(245,1);
nstot=zeros(245,1);
for n=14:-1:2
    igood=find(~isnan(SBvol(n,:)));
    if n > 2 & n < 9
        634+n;
        stot(igood)=stot(igood)+SBvol(n,igood)';
        nstot(igood)=nstot(igood)+1;
    end
    if n == 2
        plot(datetime(SBvolDates(igood),'convertfrom','datenum'),SBvol(n,igood),...
        '*-','color','m','linewidth',2,'displayname',['Mop ' num2str(634+n)])
    else
        plot(datetime(SBvolDates(igood),'convertfrom','datenum'),SBvol(n,igood),...
        '*-','color',col(n-2,:),'linewidth',2,'displayname',['Mop ' num2str(634+n)])
        hold on;
    end
end

igood=find(stot > 0 & nstot == max(nstot));
plot(datetime(SBvolDates(igood),'convertfrom','datenum'),stot(igood),...
        '*:','color','k','linewidth',2,'displayname','\Sigma637-42')


grid on;
xtickformat('MMM-dd');set(gca,'fontsize',16,'color',[.8 .8 .8])
ylabel('MSL Beach Vol (m^{3}/m-shoreline)')
lg=legend('location','eastoutside');lg.Color=([ .8 .8 .8]);
set(gca,'xlim',[datetime(2023,1,1) datetime(2023,8,1)],'ylim',[0 500])
title('2023 South End of Pad, MOP 637 = Inlet')
set(gcf,'InvertHardcopy','off')
makepng('SolanaSouthEndVolumesSxyMSL2023.png')
%%
figure('position',[95          55        1090         742]);
col=flipud(jet(12));

ax1=axes('position',[0.075 0.05 .9 .4]);
for n=31:-1:21
        stn=['D' num2str(634+n,'%4.4i')];
        [wavetime,Sxy]=GetMopSxy(stn);
   
            plot(wavetime,Sxy,...
        '-','color',col(n-19,:),'linewidth',1,'displayname',num2str(634+n))
        hold on;
end
grid on;
set(gca,'xlim',[datetime(2023,1,1) datetime(2023,8,1)])
xtickformat('MMM-dd');set(gca,'fontsize',16,'color',[.8 .8 .8])
ylabel('Tot Sxy (m^2) + = North')
lg=legend('location','eastoutside');lg.Color=([ .8 .8 .8]);

ax2=axes('position',[0.075 0.55 .9 .4]);

for n=31:-1:21
    igood=find(~isnan(SBvol(n,:)));
    % if n == 2
    %     plot(datetime(SBvolDates(igood),'convertfrom','datenum'),SBvol(n,igood),...
    %     '*-','color','m','linewidth',2,'displayname',num2str(634+n))
    % else
        plot(datetime(SBvolDates(igood),'convertfrom','datenum'),SBvol(n,igood),...
        '*-','color',col(n-19,:),'linewidth',2,'displayname',num2str(634+n))
        hold on;
    %end
end
grid on;
xtickformat('MMM-dd');set(gca,'fontsize',16,'color',[.8 .8 .8])
ylabel('MSL Beach Volume (m^{3}/m-shoreline)')
lg=legend('location','eastoutside');lg.Color=([ .8 .8 .8]);
set(gca,'xlim',[datetime(2023,1,1) datetime(2023,8,1)],'ylim',[0 300]);
title('North End of Pad : Mop 665 = Tabletop Reef')
set(gcf,'InvertHardcopy','off')
makepng('SolanaNorthEndVolumesSxyMSL2023.png')

% figure;plot(X1Dcpg,TZ1Dglo,'.-')
% hold on;grid on;
% plot(X1Dcpg,TZ1Dmon(end,:),'r.-')
% plot(X1Dcpg,TZ1Dqtr(end-1,:),'g.-')
% plot(X1Dcpg,TZ1Dsea(end,:),'c.-')
% plot(X1Dcpg,TZ1Dann(end,:),'k.-')
% 
% figure;hold on;grid on
% for n=1:numel(TZdatetime.mon)
%     plot3(TZdatetime.mon(n)+X1Dcpg*0,X1Dcpg,TZ1Dmon(n,:),'g-','linewidth',2)
% end
% 
% %figure;hold on;grid on
% for n=1:numel(TZdatetime.qtr)
%     plot3(TZdatetime.qtr(n)+X1Dcpg*0,X1Dcpg,TZ1Dqtr(n,:),'r-','linewidth',2)
% end
% 
% 
% %figure;hold on;grid on
% for n=1:numel(TZdatetime.sea)
%     plot3(TZdatetime.sea(n)+X1Dcpg*0,X1Dcpg,TZ1Dsea(n,:),'c-')
% end
% 
% 
% %figure;hold on;grid on
% for n=1:numel(TZdatetime.ann)
%     plot3(TZdatetime.ann(n)+X1Dcpg*0,X1Dcpg,TZ1Dann(n,:),'k-','linewidth',2)
% end
