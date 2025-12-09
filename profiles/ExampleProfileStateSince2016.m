% Example profile analysis focused on the the state of a Mop profile at
% time A relative to its monthly long-tmer mean and min,max,std

close all
clearvars

NumSubTrans=11;
XgapTol=15;
YdistTol=25;

nm=0;
for m=530:780
    m
 nm=nm+1;
 
 load(['M' num2str(m,'%5.5i') 'SA.mat'],'SA');
sy=year(datetime([SA.Datenum],'ConvertFrom','datenum'));
sm=month(datetime([SA.Datenum],'ConvertFrom','datenum'));

%sdx=find([SA.Datenum] >= datenum(2015,1,1));
sdx=find( sm == 8 );% datenum(2015,1,1));
%sdx=find( sy == 2020 | sy == 2024 );
SA=SA(sdx);

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
Zdatetime(end)

%figure;
%imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmean,'alphadata',~isnan(Z1Dtrans));
%datetick('y');
% pcolor(X1Dcpg,Zdatetime,Z1Dmean);shading flat;
% BeachColorbar;

[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
    TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%

Znavd88=1.344;
[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,TZ1Dmon); 

MopNum(nm)=m;
Zt(nm).mon=TZdatetime.mon;
BWmin(nm).mon=BeachWidthMin;
BWmax(nm).mon=BeachWidthMax;

% figure;plot(TZdatetime.mon,BeachWidthMin,'+-',TZdatetime.mon,BeachWidthMax,'o-')
% title('Monthly Mean Beach Width');grid on;

XbackgapTol=5;
[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,TZ1Dmon);
BVmin(nm).mon=BeachVolMin;
BVmax(nm).mon=BeachVolMax;

end

save ProfileWidthsVolumes.mat Zt BWmin BWmax BVmin BVmax
%%
figure;plot(Zt(100).mon,BVmin(100).mon,'+-')
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
