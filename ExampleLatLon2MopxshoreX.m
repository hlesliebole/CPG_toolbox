
% Example code to convert Lat,Lon point to the closest Mop (Nmop) subtransect
%  xshore profile distance (X) value(rounded to the nearest meter).

% 1. find nearest mop area to the Lat, Lon point.  
% 2. projects point onto the nearest subtransect in the mop area
%    to calculate its X distance from the Mop back beach line
%    (line connecting Mop back beach points).  
%
%  Number of subtransects is set to 101 in code (~1m alongshore spacing)
%  and returns fractional mop numbers to 2 decimal places.
%  

% set path to CPG MOP files from your machine
addpath /Volumes/group/MOPS  % folder with MOP mat files
addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% set example location in TP nearshore waters
Lat=32.93;
Lon=-117.261;
% Lon=-117.1220;
% Lat=32.5339;

% convert to UTM
[Xutm,Yutm,UTMzone]=deg2utm(Lat,Lon);

% find nearest mop area to utm point
Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);

% load Mop point info table
load MopTableUTM.mat

% divide mop area into 101 mop subtransects (~1m apart alongshore) 
% with 1m xshore resolution, with an extra 100m of back beach for each
% in case the point falls landward of the back beach line
Ntransects=101;
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,Nmop,Ntransects,[-100 0]);

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
  [row,col] = ind2sub(size(xst),NearIdx); 

% fractional mop number based on nearest subtransect number
FracMopNumber=Nmop+(row-(Ntransects+1)/2)/(Ntransects-1);
% xshore distance along the subtransect
X=x1d(col);

fprintf('Input Lat,Lon: %10.6f %10.6f\n',Lat,Lon)
fprintf('Fractional Mop Number: %10.3f\n',FracMopNumber)
fprintf('Xshore Distance from Back Beach Line (m): %i\n',X)

