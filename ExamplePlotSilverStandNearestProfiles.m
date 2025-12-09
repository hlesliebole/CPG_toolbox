% addpath /Volumes/mops/toolbox
% addpath /Volumes/mops

figure('position',[204 99 1120 637]);
ax1=axes;hold on;
col=jet(numel(105:117));

n=0;nn=0;
for MopNumber=105:117
n=n+1;
matfile=sprintf('M%5.5iSA.mat',MopNumber);

eval(['load ' matfile]);
Ytol=50;
Xtol=10;
%[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);
[X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanNearestProfiles(SA);

%% show trk back beach point on map
%[Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X0BeachOnly);
%plot(Lon,Lat,'m*','markersize',20)

idx=find([SA.Datenum] == datenum(2010,2,26));
if numel(idx) > 1;idx=idx(1);end
nn=nn+1;
pl(nn)=plot(X1Dt-X0BeachOnly,Zf(idx,:),'-','color',col(n,:),'linewidth',2,...
    'DisplayName',['Mop ' num2str(MopNumber) ' 26 Feb 2010']);
nn=nn+1;
pl(nn)=plot(X1Dt-X0BeachOnly,Zf(end,:),'--','color',col(n,:),'linewidth',1,...
    'DisplayName',['Mop ' num2str(MopNumber) ' 21 Jul 2023']);

end
set(gca,'fontsize',14,'xdir','reverse','xlim',[0 40],'color',[.7 .7 .7]);grid on;
legend(pl,'location','eastoutside');
xlabel('Distance from Parking Lot Edge (m)')
ylabel('Elev (m, Navd88)')
set(gcf,'inverthardcopy','off')
makepng('SilverStrandProfiles.png');

%plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)

