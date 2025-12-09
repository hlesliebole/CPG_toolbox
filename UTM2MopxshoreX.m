function [Nmop,X]=UTM2MopxshoreX(Xutm,Yutm)
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

%[Xutm,Yutm,UTMzone]=deg2utm(Lat,Lon);

% find nearest mop area
Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);

% load Mop point info table
load MopTableUTM.mat

for mop=unique(Nmop) % loop though range of mop numbers
    
    n=find(Nmop == mop); % points closest to this mop number
    
% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach and 500m offshore for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,mop,20,[-100 500]);

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm(n)),double(Xutm(n))],'euclidean','smallest',1);

% define X based on the nearest transect line point
[row,col] = ind2sub(size(xst),NearIdx);
X(n)=x1d(col);

end

end