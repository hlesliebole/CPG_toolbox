% Example profile analysis focused on the the state of a Mop profile at
% time A relative to its monthly long-tmer mean and min,max,std

close all
clearvars

NumSubTrans=11;
XgapTol=15;
YdistTol=25;
XbackgapTol=10;

Mop1=576;Mop1Name='MOP 576';
cs1='usa_CA_0011-0008';
%Mop1=506;Mop1Name='Mop 506: The Hole';
%cs1='usa_CA_0009-0010';
% Mop1=498;Mop1Name='Mop 498: Beach & Tennis Club';
% cs1='usa_CA_0009-0002';
% Mop1=496;Mop1Name='Mop 496: Marine Room';
% cs1='usa_CA_0009-0000';
Mop2=584;Mop2Name='MOP 584'
cs2='usa_CA_0009-0016';
%Mop2=511;Mop2Name='Mop 511: SIO Lifeguard Station';
%cs2='usa_CA_0009-0015';
matfile1=['M' num2str(Mop1,'%5.5i') 'SA.mat'];
matfile2=['M' num2str(Mop2,'%5.5i') 'SA.mat'];

load(matfile1)
[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

% figure;
% %imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmean,'alphadata',~isnan(Z1Dtrans));
% %datetick('y');
% pcolor(X1Dcpg,Zdatetime,Z1Dmean);shading flat;
% BeachColorbar;

[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
    TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%

Znavd88=1.344;
[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,TZ1Dmon); 

[VosDatetimes,VosX]=GetCoastSat(cs1);
% CoastSat anomaly relative to survey mean MHW
VosX=(VosX-mean(VosX))+mean(BeachWidthMin,'omitnan');


f1=figure('position',[ 1          26        1440         771]);
subplot(2,1,1);
p0=plot(VosDatetimes,VosX,'.-','color',[.6 .6 1]);hold on;
p1=plot(TZdatetime.mon,BeachWidthMin,'b+-','linewidth',2);
plot(TZdatetime.mon,BeachWidthMax,'bo-','linewidth',2);
title('Monthly Mean MHW Beach Widths (Thick = Survey ; Thin = CoastSat)');grid on;hold on;
% set(gca,'ylim',[-5 65],'xlim',[datetime(1984,1,1) datetime(2025,1,1)],...
%     'xtick',datetime(1984:1:2026,1,1),'fontsize',16)
set(gca,'ylim',[-5 Inf],'xlim',[datetime(1984,1,1) datetime(2025,1,1)],...
    'xtick',datetime(1984:1:2026,1,1),'fontsize',16)
ylabel('MHW Width (m)')



[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,TZ1Dmon);
%f2=figure;
subplot(2,1,2);
plot(TZdatetime.mon,BeachVolMin,'b+-','linewidth',2);
plot(TZdatetime.mon,BeachVolMin,'bo-','linewidth',2)
title('Monthly Mean MHW Beach Volumes');grid on;hold on;
set(gca,'xlim',[datetime(1984,1,1) datetime(2025,1,1)],...
    'xtick',datetime(1984:1:2026,1,1),'fontsize',16)
ylabel('MHW Vol (m^{3}/m-shore)')


%% -----
load(matfile2)


idx=find([SA.Datenum]==datenum(2024,1,16));SA(idx)=[];
idx=find([SA.Datenum]==datenum(2024,1,24));SA(idx)=[];
idx=find([SA.Datenum]==datenum(2024,1,30));SA(idx)=[];
idx=find([SA.Datenum]==datenum(2024,2,12));SA(idx)=[];

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
    TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
Znavd88=1.344;
[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,TZ1Dmon); 
[BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dcpg,TZ1Dmon);

[VosDatetimes,VosX]=GetCoastSat(cs2);
% CoastSat anomaly relative to survey mean MHW
VosX=(VosX-mean(VosX))+mean(BeachWidthMin,'omitnan');

subplot(2,1,1);
p0=plot(VosDatetimes,VosX,'.-','color',[1 .6 .6]);hold on;
p2=plot(TZdatetime.mon,BeachWidthMin,'r+-','linewidth',2);
plot(TZdatetime.mon,BeachWidthMax,'ro-','linewidth',2);
legend([p1(1) p2(1)],Mop1Name,Mop2Name,'location','northwest')
subplot(2,1,2);
plot(TZdatetime.mon,BeachVolMin,'r+-','linewidth',2);
plot(TZdatetime.mon,BeachVolMin,'ro-','linewidth',2)

function [VosDatetimeMon,VosXmon]=GetCoastSat(TransectName)

%TransectName=['usa_CA_0011-0016'];

 % get CoastSat transect shoreline time series from website
 url=['http://coastsat.wrl.unsw.edu.au/time-series/' TransectName '/'];
  s=webread(url);

ss=strsplit(s,'\n');

% add any valid data to time series
k=0; % data counter
for n=1:size(ss,2)
    if ~isempty(ss{n})
       st=strsplit(regexprep(ss{n},',',' '),' ');
       k=k+1;
       VosDatetimes(k)=datetime([st{1} ' ' st{2}]);
       VosX(k)=str2double(st{3});
    end
end

% check for and remove any time series shoreline position NaNs
idx=find(isnan(VosX));
if numel(idx) > 0
    fprintf(' *** Removing %i CoastSat shoreline "None"s\n',numel(idx));
    VosDatetimes(idx)=[];
    VosX(idx)=[];
end

fprintf('CoastSat Valid N= %i\n',numel(VosX));

% survey datetimes
Zdt=VosDatetimes;

%% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) > 0
           n=n+1;
           Zymdate(n)=datetime(y,m,15,0,0,0);
           
           if numel(idx) == 1
               VosXmon(n)=VosX(idx);
           elseif numel(idx) > 1
               VosXmon(n)=mean(VosX(idx),'omitnan');
           end
       end
   end
end

VosDatetimeMon=Zymdate;

end
