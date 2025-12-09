% Example profile analysis focused on the the state of a Mop profile at
% time A relative to its monthly long-tmer mean and min,max,std

close all
clearvars

NumSubTrans=11;
XgapTol=5;
YdistTol=2;
XbackgapTol=15;
Znavd88=1.344;%0.774;%

LastN=2; % keep the last 2

MopStart=567;%635;%567%2
MopEnd=636;%684;%636%926
nmops=numel(MopStart:MopEnd);
nm=0;
% BWmin=NaN(1,nmops);
% BWmax=NaN(1,nmops);
% BVmin=NaN(1,nmops);
% BVmax=NaN(1,nmops);
MopNum=MopStart:MopEnd;
for m=MopStart:MopEnd
%    for m=530:780
    fprintf('%i of %i\n',m,MopEnd)
 nm=nm+1;
 Zt(nm,1:LastN)=NaN;
 BWmin(nm,1:LastN)=NaN;
 BWmax(nm,1:LastN)=NaN;
 BVmin(nm,1:LastN)=NaN;
 BVmax(nm,1:LastN)=NaN;
 
load(['M' num2str(m,'%5.5i') 'SA.mat'],'SA');

ldx=find(strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR'));

if numel(ldx) >= LastN

sy=year(datetime([SA.Datenum],'ConvertFrom','datenum'));
sm=month(datetime([SA.Datenum],'ConvertFrom','datenum'));

%sdx=find([SA.Datenum] >= datenum(2015,1,1));
%sdx=find( sm == 8 );% datenum(2015,1,1));
%sdx=find( sy == 2020 | sy == 2024 );
% most recent 2
SA=SA(ldx(end-10:end));

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

%figure;
%imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmean,'alphadata',~isnan(Z1Dtrans));
%datetick('y');
% pcolor(X1Dcpg,Zdatetime,Z1Dmean);shading flat;
% BeachColorbar;
% 
% [TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
%     TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%


[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,Z1Dmedian); 


% figure;plot(TZdatetime.mon,BeachWidthMin,'+-',TZdatetime.mon,BeachWidthMax,'o-')
% title('Monthly Mean Beach Width');grid on;


[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,Z1Dmedian);

MopNum(nm)=m;
for nls=1:LastN
    %ns=numel(ldx)-LastN+nls;
    ns=11-LastN+nls;
 Zt(nm,nls)=datenum(Zdatetime(ns));
 BWmin(nm,nls)=BeachWidthMin(ns);
 BWmax(nm,nls)=BeachWidthMax(ns);
 BVmin(nm,nls)=BeachVolMin(ns);
 BVmax(nm,nls)=BeachVolMax(ns);
end

end
end

save LidarLast2ProfileWidthsVolumes.mat MopNum Zt BWmin BWmax BVmin BVmax
%%
% btimes=vertcat(Zt.last2);
% bwmin=vertcat(BWmin.last2);
% bwmax=vertcat(BWmax.last2);
% bvmin=vertcat(BVmin.last2);
% bvmax=vertcat(BVmax.last2);
% 
% figure;plot(MopNum,bwmin(:,1),'+-',MopNum,bwmin(:,2),'o-')
% %figure;plot(btimes(:,1),bwmin(:,1),'+-',btimes(:,2),bwmin(:,2),'o-')
% title(' MHW Beach Vol');grid on;
% 

% unuque most recent survey dates
udates=unique(Zt(~isnan(Zt(:,LastN)),LastN));udates(end)=[];
%datestr(unique(Zt(~isnan(Zt(:,2)),2)))

figure('position',[23  218  1350   700]);
subplot(2,1,1)
%dBW=BWmin(:,2)-BWmin(:,1);
idx=find(Zt(:,2) == udates(end));
bar(MopNum(idx),[BWmin(idx,1) BWmin(idx,2)]);grid on;hold on;
legend(datestr(Zt(idx(1),1)),datestr(Zt(idx(1),2)))
xl=get(gca,'xlim');
set(gca,'xtick',round(xl(1)):round(xl(2)),'fontsize',16);
ylabel('Beach Width (m)')
title(['MOP Beach Width (m) @ z= ' num2str(Znavd88,'%5.3f') ' NAVD88'] )
%bar(MopNum(idx),BWmin(idx,2));grid on;hold on;
% idx=find(dBW < 0 & Zt(:,2) == udates(end));
% bar(MopNum(idx),dBW(idx),'r');
subplot(2,1,2)
dBW=BWmin(:,2)-BWmin(:,1);
idx=find(dBW > 0 & Zt(:,2) == udates(end));
bar(MopNum(idx),dBW(idx),'FaceColor',[0 .7  0 ]);grid on;hold on;
idx=find(dBW < 0 & Zt(:,2) == udates(end));
bar(MopNum(idx),dBW(idx),'r');
set(gca,'xlim',xl);
set(gca,'xtick',round(xl(1)):round(xl(2)),'fontsize',16);
ylabel('Beach Width Change(m)')


%bar(MopNum(idx),BVmin(idx,2),'FaceColor',[0 .7  0 ]);grid on;hold on;
%%
figure('position',[23  218  1350   700]);
subplot(2,1,1)
%dBW=BWmin(:,2)-BWmin(:,1);
idx=find(Zt(:,2) == udates(end));
bar(MopNum(idx),[BVmin(idx,1) BVmin(idx,2)]);grid on;hold on;
legend(datestr(Zt(idx(1),1)),datestr(Zt(idx(1),2)))
xl=get(gca,'xlim');
set(gca,'xtick',round(xl(1)):round(xl(2)),'fontsize',16);
ylabel('Beach Volume (m^{3}/m)')
title(['MOP Beach Volumes (m^{3}/m) above z= ' num2str(Znavd88,'%5.3f') ' NAVD88'] )
%bar(MopNum(idx),BWmin(idx,2));grid on;hold on;
% idx=find(dBW < 0 & Zt(:,2) == udates(end));
% bar(MopNum(idx),dBW(idx),'r');
subplot(2,1,2)
dBV=BVmin(:,2)-BVmin(:,1);
idx=find(dBV > 0 & Zt(:,2) == udates(end));
bar(MopNum(idx),dBV(idx),'FaceColor',[0 .7  0 ]);grid on;hold on;
idx=find(dBV < 0 & Zt(:,2) == udates(end));
bar(MopNum(idx),dBV(idx),'r');
set(gca,'xlim',xl);
set(gca,'xtick',round(xl(1)):round(xl(2)),'fontsize',16);
ylabel('Beach Volume Change(m^{3}/m)')



% figure('position',[23         218        1350         466]);
% plot(MopNum,BVmin(:,2)-BVmin(:,1),'-');grid on
% figure('position',[23         218        1350         466]);
% plot(MopNum,bvmin(:,2)-bvmin(:,1),'-');grid on


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
