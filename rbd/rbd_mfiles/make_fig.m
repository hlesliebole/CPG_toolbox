%
% Makes a color dot plot of nos bathymetry data read in by
% load_rbd.m
%
% mfile to make the figure window for the bathymetry plot
%  works on 15 inch Sun monitor

% Scale the axes appropriately

%lon_min=floor(min(lon));
%lat_min=floor(min(lat));
%lon_max=ceil(max(lon));
%lat_max=ceil(max(lat));
lon_min=(min(lon));
lat_min=(min(lat));
lon_max=(max(lon));
lat_max=(max(lat));

avlat=(lat_max+lat_min)/2;
avlat=avlat*(pi/180);

% aspect factor for lon vs. lat  = cos(deg lat)
fac=cos(avlat);

% specify size of the figure window

dx=(lon_max-lon_min)*fac;
dy=lat_max-lat_min;
dfx=.05*dx;
dfy=.05*dy;

% adjust rect below if figure doesn't fit on screen

L=750;
if dx > dy
top=200+L*dy/dx;
rect=[275 30 L top];
center=[ L/2 top/2];
else
width=L;
top=200+L;
rect=[275 30 width top];
center=[ width/2 top/2 ];
end

dy/dx
rect
h=figure('Numbertitle','off','Position',rect,...
          'name','Bathymetry Editing Program');

rat=[dx/dy fac];
rect=[0.12 0.10 .8 .8];
limy=[lat_min-dfy lat_max+dfy];
limx=[lon_min-dfx lon_max+dfx];
h=axes('Position',rect,'Box','on','DataAspectRatio',...
        [rat 1],'Xlim',limx,'Ylim',limy);
title(plot_ttl);
xlabel('Longitude')
ylabel('Latitude')
grid on
hold;
