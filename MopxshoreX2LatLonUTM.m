function [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNumber,X)

% Converts Mop xshore profile X to Lat Lon and UTM coords

% Calls utm2deg.m


% NOTE: hardwired for UTM zome '11 S'
UTMzone='11 S';

% load Mop point info table
load MopTableUTM.mat

% first convert xshore x value to utm coordinates

% mop transect angle in utm coord frame
MopTransectAngle=atan((Mop.BackYutm(MopNumber)-Mop.OffYutm(MopNumber))...
    ./(Mop.BackXutm(MopNumber)-Mop.OffXutm(MopNumber)));

% convert X to xutm, yutm
Xutm=Mop.BackXutm(MopNumber)-X*cos(MopTransectAngle);
Yutm=Mop.BackYutm(MopNumber)-X*sin(MopTransectAngle);

% convert to lat lon using utm2deg.m
[Lat,Lon]=utm2deg(Xutm,Yutm,repmat(UTMzone,[length(Xutm) 1]));

end