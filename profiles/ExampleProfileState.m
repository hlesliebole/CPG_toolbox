% Example profile analysis focused on the the state of a Mop profile at
% time A relative to its monthly long-tmer mean and min,max,std

% close all
% clearvars
%% Configuration
mpath = '/volumes/group/MOPS/';  % Path to CPGMOP data files
addpath('/volumes/group/MOPS');
addpath('/volumes/group/MOPS/toolbox/');

NumSubTrans=51;
XgapTol=15;
YdistTol=25;

load M00649SA.mat

%SA=SA(end-1:end);

% idx=find([SA.Datenum]==datenum(2024,1,16));SA(idx)=[];
% idx=find([SA.Datenum]==datenum(2024,1,24));SA(idx)=[];
% idx=find([SA.Datenum]==datenum(2024,1,30));SA(idx)=[];
% idx=find([SA.Datenum]==datenum(2024,2,12));SA(idx)=[];

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

figure;
%imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmean,'alphadata',~isnan(Z1Dtrans));
%datetick('y');
pcolor(X1Dcpg,Zdatetime,Z1Dmean);shading flat;
BeachColorbar;

[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
    TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%

Znavd88=1.344;
[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,TZ1Dmon); 

figure;plot(TZdatetime.mon,BeachWidthMin,'+-',TZdatetime.mon,BeachWidthMax,'o-')
title('Monthly Mean Beach Width');grid on;

XbackgapTol=5;
[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,TZ1Dmon);
figure;plot(TZdatetime.mon,BeachVolMin,'+-',TZdatetime.mon,BeachVolMin,'o-')
title('Monthly Mean MHW Beach Vol');grid on;


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
