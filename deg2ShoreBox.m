function [Xsb,Ysb]=deg2ShoreBox(Lat,Lon,Mop1,Mop2)

%--------------------------------------------------------------
% Dependencies:
%
%   Loads the CA Mop transect information in matlab table variable "Mop"
%      from MopTable.mat
%
%--------------------------------------------------------------

% input:
% 
% Lat = Latitudes
% Lon = Longitudes
% Mop1 = downcoast MOP number or station name 
% Mop2 = upcoast MOP number or station name  

% Mop1 and Mop2 can be numeric mop numbers (CA all numeric convention 
%  starting at US-MX border) or the character string CDIP MOP 
%  identifier.  eg. 663 or ?D0663?

% output:
%  
% Xsb = shorebox x values
% Ysb = shorebox y values
% 

% if Mop numbers are entered as character string names, get numeric
% equivalent
if(ischar(Mop1))
    Mop1=find(cellfun('isempty',strfind(Mop.Name,Mop1)) == 0);
    Mop2=find(cellfun('isempty',strfind(Mop.Name,Mop2)) == 0);
end

% load mop table
load MopTable.mat Mop

% reduce table to desired mop range
MopSB=Mop(Mop1:Mop2,:);

% convert mop transect back beach and offshore points to utm
[MopSB.BackLon,MopSB.BackLat,utmzone]=deg2utm(MopSB.BackLat,MopSB.BackLon);
[MopSB.OffLon,MopSB.OffLat,utmzone]=deg2utm(MopSB.OffLat,MopSB.OffLon);

% initial axis and (0,0) origin based on first and last Mop
x0=MopSB.BackLon(1);y0=MopSB.BackLat(1);
xe=MopSB.BackLon(end);ye=MopSB.BackLat(end);
theta0=atan2((ye-y0),(xe-x0)); % angle of axis

% transform all the Mop backbeach points
d=sqrt((MopSB.BackLon-x0).^2+(MopSB.BackLat-y0).^2);
theta=atan2((MopSB.BackLat-y0),(MopSB.BackLon-x0));
thetat=(theta-theta0); % rotation
MopSB.BackLon=d.*cos(thetat); % rotated mop backbeach x values 
MopSB.BackLat=d.*sin(thetat); % rotated mop backbeach y values

% transform all the Mop offshore points
d=sqrt((MopSB.OffLon-x0).^2+(MopSB.OffLat-y0).^2);
theta=atan2((MopSB.OffLat-y0),(MopSB.OffLon-x0));
thetat=(theta-theta0); % rotation
MopSB.OffLon=d.*cos(thetat); % rotated mop backbeach x values 
MopSB.OffLat=d.*sin(thetat); % rotated mop backbeach y values

% Now make final axis origin by shifting y axis values so mininum mop
% back beach y value = 0
dy=min(MopSB.BackLat);

% apply y shift to shorebox MopSB y coords
MopSB.BackLat=MopSB.BackLat-dy;
MopSB.OffLat=MopSB.OffLat-dy;

% now transform survey data with the x0,y0,theta0 and dy transformation
%  parameters

% convert lat lon data to utm xdu ydu
[xdu,ydu,utmzone]=deg2utm(Lat,Lon);

d=sqrt((xdu-x0).^2+(ydu-y0).^2); % distance from shorebox origin
theta=atan2((ydu-y0),(xdu-x0)); % angle from shorebox origin
thetat=(theta-theta0); % rotation angle
Xsb=d.*cos(thetat); % rotate x values 
Ysb=d.*sin(thetat); % rotate y values

Ysb=Ysb-dy; % apply y shift to y coords


end