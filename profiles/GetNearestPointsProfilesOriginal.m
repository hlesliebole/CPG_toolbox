function [X1D,Z1D]=GetAllNearestPointsProfilesOriginal(SA,dy,dx)

%  Calculates nearest points-based profiles for all the surveys in 
%  the input SA struct array file
% eg. for the median profile shape of TPines surveys
%   load M00583SA.mat
%   [X1D,Z1D]=GetNearestPointsProfile(SA,25,5);
%   figure;plot(X1D,median(rmoutliers(Z1D,"mean"),'omitnan'));grid on
%   figure;plot(median(rmoutliers(Z1D,"mean"),'omitnan'),gradient(median(rmoutliers(Z1D,"mean"),'omitnan')),'*');grid on

% Input:
%   SA = 1m spatial avg struct array data file 
%   dy = alongshore distance from mop line tolerance (m)
%   dx = max xshore dist tolerance for small gap filling interpolation (m)

%  Returns:
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1D(size(SA,2),numel(X1D)) = 2d matrix of nearest profiles for all SA surveys

% Dependencies:
%
%  needs /MOPS and MOPS/toolbox in its path
%  calls GetTransectLines.m
%  uses  MopTableUTM.mat

% Rather than interpolating a profile from gridded survey points,
% this code 
%
%  1) loads the 1m spatial averaged survey point data and
%  2) finds the distance of each from the mop line defined 
%      as a line of xshore points with 1m spacing.
%  3) for each xshore point that has one or more survey points
%      that are closer to it than the other mop line points,
%      find the closest of those survey points 
%
% The resulting profile variables are
%
%     Xuniq = the unique xshore distances with survey data (1m integers)
%     Znear = the elevation of the survey point nearest to each Xuniq point (m,
%             navd88)
%     Zdist = the distance of each Znear point from the mop line (m)
%
% It then makes mop line profile vectors x1d and z1di that use a 
%   distance from the mop line threshhold, dy, to define valid profile 
%   Znear data, and a maximum allowable spacing between Xuniq points, dx, 
%   that can be filled using 1d interpolation of the remaining Znear points,
%   if interpolated/regular spaced xshore profile info is needed.
 

% set path to CPG MOP files from your machine
% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');
 
%  Get the main transect line as UTM x,y points ~1m apart
%   This also sets the fixed x1d transect point starting
%   -100m behind the mop back beach point out to the offshore point
[X1D,xt,yt,xst,yst]=GetTransectLines(Mop,SA(1).Mopnum,1,[-100 0]);

Z1D=NaN(size(SA,2),numel(X1D));

% loop through all surveys

for nsurv=1:size(SA,2)

clear Zdist Znear

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(SA(nsurv).Y),double(SA(nsurv).X)],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xst),NearIdx); 

% xshore distance (1m xshore resolution) along the Mop transect for each survey point
X=X1D(col);

% For the xshore distances (1m xshore resolution) with data (Xuniq), 
%  find the nearest survey point elevation for each (Znear) and
%  how far away from the line it was (Zdist).

Xuniq=unique(X);

n=0;
for x=Xuniq
    n=n+1; % xshore point counter
    idx=find(X == x);
    [Zdist(n),imin]=min(dp(idx));
    Znear(n)=SA(nsurv).Z(idx(imin));   
end

% make an interpolated profile based on distance from 
% mop line (dy) and xshore spatial gap (dx) tolerances

% set nearest points with Zdist > dy to NaNs
Znear(Zdist > dy)=NaN;

% seed the regular spaced profile vectors (x1d and z1d) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(X1D,Xuniq));
z1d=X1D*NaN;z1d(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(z1d);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= dx are removed. NaN
%  no data points in the gaps larger than dx are kept.

x1dv=X1D;x1dv(isnan(z1d) & sz <= dx)=[];
z1dv=z1d;z1dv(isnan(z1d) & sz <= dx)=[];

% interpolate to fill the small gaps
%z1di=interp1(x1dv,z1dv,x1d);
Z1D(nsurv,:)=interp1(x1dv,z1dv,X1D);
end

%figure;pcolor(X1D,1:size(SA,2),Z1D);shading flat

end



