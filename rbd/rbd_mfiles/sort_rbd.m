%
% matlab script to sort the data by latitude
%

' sorting data by latitude '

[Y,I]=sort(lat);
lat=Y;
Y=lon(I);
lon=Y;
Y=dep(I);
dep=Y;
clear Y I;

