% load cf data
load SANDAG_Transect_Data.mat
% plot the locations of the first x930 transect in the struct array
Xutm=Transect_Data.x930.XUTM(:,1);
Yutm=Transect_Data.x930.YUTM(:,1);
Znavd88=Transect_Data.x930.Znavd88(:,1);

% convert the utm coords to lat on
[lat,lon]=utm2deg(Xutm,Yutm,repmat('11 S',[numel(Xutm) 1]));

plot(lon,lat,'+');hold on;

% overlay on google map
plot_google_map