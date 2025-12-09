load M00654SA.mat
NumSubTrans=11;
XgapTol=5;
YdistTol=5;

% [X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
%   GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

%%

Znavd88=1.344;
[BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,Z1Dtrans); 

figure;plot(Zdatetime,BeachWidthMin)

%%
idx=find(Zdatetime == datetime(2018,10,9),1,'first');
clear pl
figure('position',[14          68        1405         725]);
axes('position',[.05 .55 .9 .4],'xdir','reverse','linewidth',2,'fontsize',16);grid on;box on;
np=0;
for n=0:7
 if contains(SA(idx+n).Source,'Trk')
     np=np+1;
  pl(np)=plot(X1Dmop,Z1Dmedian(idx+n,:)-0.774,'linewidth',2,'DisplayName',string(Zdatetime(idx+n)));
  pdatetime(np)=Zdatetime(idx+n);
 end
hold on;
end

np=np+1;
pl(np)=plot(X1Dmop,Z1Dmedian(end,:)-0.774,'m-','linewidth',2,'DisplayName',strcat(string(Zdatetime(end)),' (for Reference)'));
xlabel('Cross-Shore Distance (m)');ylabel('Elevation (m, MSL)')
xl=get(gca,'xlim');
%plot(xl,[.774 .774],'k--','linewidth',2);text(xl(2),.9,' MSL','fontsize',16,'linewidth',2);
legend(pl,'location','northwest');
title('Erosion of Fletcher MOP 654 2018 Nourishment Pad')

addpath ../../../runup
ax2=axes('position',[.05 .075 .9 .35]);

Rstart=datenum(2018,10,1);  % runup time series start date
Rend=datenum(2019,3,1); % runup time series end date
Bslope=0.1;

PlotStockdonHindcast
hold on;
tl=plot(moptime,moptime*0+3.2,'k--','DisplayName','2018 Terrace Elev');

yl=get(gca,'ylim');
for n=1:np-1
    xt=datenum(pdatetime(n));
    plot([xt xt],yl,'-','color',pl(n).Color,'linewidth',2)
end

legend(tl)

makepng('Fletcher2018PadErosion.png')

