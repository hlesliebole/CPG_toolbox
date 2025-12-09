function [Lat,Lon]=ShoreBox2deg(Xsb,Ysb,Mop1,Mop2)

% Transforms shorebox coordinates back to Lat Lon
% input:
% 
% Xsb = shorebox x values
% Ysb = shorebox y values
% Mop1 = downcoast MOP number or station name 
% Mop2 = upcoast MOP number or station name  

% Mop1 and Mop2 can be numeric mop numbers (CA all numeric convention 
%  starting at US-MX border) or the character string CDIP MOP 
%  identifier.  eg. 663 or ?D0663?

% output:
%  
% Lat = Latitudes
% Lon = Longitudes
% 

% if Mop numbers are entered as character string names, get numeric
% equivalent
if(ischar(Mop1))
    Mop1=find(cellfun('isempty',strfind(Mop.Name,Mop1)) == 0);
    Mop2=find(cellfun('isempty',strfind(Mop.Name,Mop2)) == 0);
end

% load mop table
load MopTable.mat

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
%MopSB.BackLat=MopSB.BackLat-dy;
%MopSB.OffLat=MopSB.OffLat-dy;

% now reverse transform survey data with the x0,y0,theta0 and dy transformation
%  parameters to standard UTM coordinates
Ysb=Ysb+dy; % undo y shift relative to utm y origin
d=sqrt(Xsb.^2+Ysb.^2); % distance from shorebox origin
theta=atan2(Ysb,Xsb); % angle from utm x axis
thetat=(theta+theta0); % rotation angle BACK to UTM 
xutm=x0+d.*cos(thetat); % rotate x values and add utm x origin
yutm=y0+d.*sin(thetat); % rotate y values and add utm y origin

% convert utm x y to lat lon 
%  first figure out the utmzone based on Mop1 lat lon
%[xutm,yutm,utmzone]=deg2utm(Mop(Mop1).BackLat,Mop(Mop1).BackLon);
[Lat,Lon]=utm2deg(xutm,yutm,repmat(utmzone(1,:),[length(xutm) 1]));

end