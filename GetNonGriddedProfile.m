function [x1d,z1di]=GetNonGriddedProfile(MopNumber,nsurv)

% Example code to plot SA survey data points projected onto
%  the Mop line and color coded by distance from the mop line

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

% set the mop number of interest
% MopNumber=743; % Beacons
% MopNumber=900; % Oside

% load the SA mat file
matfile=sprintf('M%5.5iSA.mat',MopNumber);
load(matfile,'SA');

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');
 
%  Get thethe main transect line as UTM x,y points ~1m apart
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,1,[-100 0]);

% select a survey to analyize
%  Use USACE SHOALS data from 26 Oct 2009 as a complex example
% nsurv=19;%22;
% nsurv=size(SA,2);%last survey;

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(SA(nsurv).Y),double(SA(nsurv).X)],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xst),NearIdx); 

% xshore distance (1m xshore resolution) along the Mop transect for each survey point
X=x1d(col);

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

% example to make an interpolated profile based on distance from 
% mop line (dy) and xshore spatial gap (dx) tolerances
dy=25; % alongshore distance from mop line tolerance
dx=5; % xshore date gap tolerance

% set nearest points with Zdist > dy to NaNs
Znear(Zdist > dy)=NaN;

% seed the regular spaced profile vectors (x1d and z1d) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(x1d,Xuniq));
z1d=x1d*NaN;z1d(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(z1d);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= dx are removed. NaN
%  no data points in the gaps larger than dx are kept.

x1dv=x1d;x1dv(isnan(z1d) & sz <= dx)=[];
z1dv=z1d;z1dv(isnan(z1d) & sz <= dx)=[];

% interpolate to fill the small gaps
z1di=interp1(x1dv,z1dv,x1d);

% smooth if mop num = 920 and n = 86
if MopNumber == 920 && nsurv == 86
    z1di = movmean(z1di, 3);
end

end



