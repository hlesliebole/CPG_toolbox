
%fl='/users/william/desktop/Jetski20230501.txt';
%fl2='/users/william/desktop/Jetski20230508.txt';
ds=[20230417
20230421
20230501
20230505
20230508
20230512];
ms=[578
    578
    578
    578
    578
    578];
me=[595
    595
    595
    595
    595
    598];

figure;hold on;

for n=1:6
fl=['/volumes/group/topobathy/' num2str(ds(n)) '_00' num2str(ms(n))...
    '_00' num2str(me(n)) '_Torrey_jetski/filtered' num2str(ds(n)) 'raw.llnezts'];
%fl2='/volumes/group/topobathy/20230508_00578_00595_Torrey_jetski/filtered20230508raw.llnezts';

%yPierDeck=[32.8669   32.8675];
yPierDeck=[32.8669   32.8671];
%yPierDeck=[32.8669   32.867075];
xPierDeck=[-117.2574 -117.2571];
% yPierDeck=[0 Inf];
% xPierDeck=[-Inf 0];
d=load(fl);
% lat=d(:,7);
% lon=d(:,8);
% z=d(:,13);
lat=d(:,1);
lon=d(:,2);
z=d(:,5);

idx=find(lat > yPierDeck(1) & lat < yPierDeck(2) &...
    lon > xPierDeck(1) & lon < xPierDeck(2) );
% 
% d=load(fl2);
% % lat2=d(:,7);
% % lon2=d(:,8);
% % z2=d(:,13);
% lat2=d(:,1);
% lon2=d(:,2);
% z2=d(:,5);

% idx2=find(lat2 > yPierDeck(1) & lat2 < yPierDeck(2) &...
%     lon2 > xPierDeck(1) & lon2 < xPierDeck(2) );

% figure;plot(lon(idx),lat(idx),'c.');hold on;plot(lon2(idx2),lat2(idx2),'m.')
% plot_google_map('MapType', 'satellite')
zmean=mean(z(idx(2:10)));
plot(z(idx(2:10)),'o','linewidth',2,'DisplayName',[num2str(ds(n)) ' : ' num2str(zmean,'%5.2f')]);hold on;
%plot(z2(idx2),'ro','linewidth',2,'DisplayName','8 May 2023');

end

set(gca,'ylim',[-7 -5],'fontsize',12);grid on;box on;
xlabel('Sample Number');ylabel('Elevation (m,NAVD88)');title('SIO Pier Launch Depths','fontsize',16);
lg=legend;
title(lg,'  Date : Mean Depth')
makepng('PierLaunchDepths.png')