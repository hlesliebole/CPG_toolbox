
%fl='/users/william/desktop/Jetski20230501.txt';
%fl2='/users/william/desktop/Jetski20230508.txt';

fl='/volumes/group/topobathy/20151216_01605_01642_newport_jumbo/filtered_clean20151216.llzts.navd88';
%fl2='/volumes/group/topobathy/20230515_00578_00598_Torrey_jetski/transects/filtered_jetski.llnezts';

%yPierDeck=[32.8669   32.8675];
yPierDeck=[32.8669   32.8671];
%yPierDeck=[32.8669   32.867075];
xPierDeck=[-117.2574 -117.2571];
yPierDeck=[0 Inf];
xPierDeck=[-Inf 0];
d=load(fl);
% lat=d(:,7);
% lon=d(:,8);
% z=d(:,13);
lat=d(:,1);
lon=d(:,2);
z=d(:,5);

idx=find(lat > yPierDeck(1) & lat < yPierDeck(2) &...
    lon > xPierDeck(1) & lon < xPierDeck(2) );

% d=load(fl2);
% % lat2=d(:,7);
% % lon2=d(:,8);
% % z2=d(:,13);
% lat2=d(:,1);
% lon2=d(:,2);
% z2=d(:,5);

% idx2=find(lat2 > yPierDeck(1) & lat2 < yPierDeck(2) &...
%     lon2 > xPierDeck(1) & lon2 < xPierDeck(2) );

figure;plot(lon(idx),lat(idx),'y.');hold on;%plot(lon2(idx2),lat2(idx2),'m.')
plot_google_map('MapType', 'satellite')

%figure;plot(z(idx),'b+');hold on;plot(z2(idx2),'m+');