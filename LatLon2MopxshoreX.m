function [Nmop,X]=LatLon2MopxshoreX(Lat,Lon)
% Converts Lat,Lon point to the closest Mop (Nmop) xshore profile 
% distance (X) value(rounded to the nearest meter in X).

% 1. find nearest mop area to the Lat, Lon point.  
% 2. projects point onto the nearest subtransect in the mop area
%    to calculate its X distance from the Mop back beach line
%    (line connecting Mop back beach point).
%
%  This is what the CPG Mop code does that creates the xshore profile
%   data in the SM struct arrays.

% Note: MOPS/toolbox needs to be in your path. Calls 
%  deg2utm.m
%  FindNearestMopTransectsUTM.m
%  GetTransectLines.m

[Xutm,Yutm,UTMzone]=deg2utm(Lat,Lon);

% find nearest mop area
Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);

% load Mop point info table
load MopTableUTM.mat

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,Nmop,20,[-100 0]);

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

% define X based on the nearest transect line point
[row,col] = ind2sub(size(xst),NearIdx);
X=x1d(col);

end