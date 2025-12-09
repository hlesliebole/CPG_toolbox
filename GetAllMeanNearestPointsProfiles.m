function [X1D,Z1Dmean,Z1Dmed,Z1Dstd]=GetAllMeanNearestPointsProfiles(SA,Ytol,Xtol,Nsubtrans)

%%

%  For each survey in the SA struct array, calculates mean and median 
%  profiles, and xshore standard deviations of a series of subtransect 
%  profiles spanning a Mop's area. 
%
%  The Mop area is divided into Nsubtrans subtransects centered on the
%  main Mop transect line. For each survey, nearest points-based profiles 
%  are calculated for each of the subtransects, and the mean,
%  median and standard deviation of the Mop area's subtransect xshore 
%  elevations are derived.
%

% Input:
%   SA = 1m spatial avg struct array data file of size(SA,2) surveys
%   Ytol = alongshore distance from mop subtransect tolerance (m)
%   Xtol = max subtransect xshore dist tolerance for small gap filling interpolation (m)
%   Nsubtrans = number of suntransects to divide the Mop area into

%  Returns:
%   X1D = vector of 1m xshore grid points (fixed range to encompass all the SA surveys)
%   Z1Dmean(size(SA,2),numel(X1D)) = 2d matrix of the mean subtransect nearest points
%                                     profiles 
%   Z1Dmed(size(SA,2),numel(X1D)) = 2d matrix of the median subtransect nearest points
%                                     profiles 
%   Z1Dstd(size(SA,2),numel(X1D)) = 2d matrix of the standard deviation of the subtransect 
%                                     nearest points profiles 

% Dependencies:
%
%  needs /MOPS and MOPS/toolbox in its path
%  calls GetTransectLines.m
%  uses  MopTableUTM.mat

%
% The nearest points profile variables are
%
%     Xuniq = the unique xshore distances with survey data (1m integers)
%     Znear = the elevation of the survey point nearest to each Xuniq point (m,
%             navd88)
%     Zdist = the distance of each Znear point from the mop line (m)
%
%   Generates makes mop line profile vectors x1d and z1di that use a 
%   distance from the mop line threshhold, Ytol, to define valid profile 
%   Znear data, and a maximum allowable spacing between Xuniq points, Xtol, 
%   that can be filled using 1d interpolation of the remaining Znear points,

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');
 
%   Get the main transect as UTM xt,yt points and subtransect lines as 
%   UTM xst,yst points ~1m apart.  This also sets the fixed 1m spaced  x1d transect 
%   points starting -100m behind the mop back beach point out to the offshore point
[X1D,xt,yt,xst,yst]=GetTransectLines(Mop,SA(1).Mopnum,Nsubtrans,[-100 0]);

Z1Dmean=NaN(size(SA,2),numel(X1D));
Z1Dmedian=NaN(size(SA,2),numel(X1D));
Z1Dstd=NaN(size(SA,2),numel(X1D));

% loop through all surveys

for nsurv=1:size(SA,2)

Z1Dsubtran=NaN(Nsubtrans,numel(X1D));

% loop through all the subtransects
for ntran=1:Nsubtrans
clear Zdist Znear

% find the distances of all dat apoints to the subtransect line points  
[dp,NearIdx]=...
    pdist2([yst(ntran,:)',xst(ntran,:)'],[double(SA(nsurv).Y),double(SA(nsurv).X)],...
    'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xst(ntran,:)),NearIdx); 

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
% mop line (Ytol) and xshore spatial gap (Xtol) tolerances

% set nearest points with Zdist > Ytol to NaNs
Znear(Zdist > Ytol)=NaN;

% seed the regular spaced profile vectors (x1d and z1d) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(X1D,Xuniq));
z1d=X1D*NaN;z1d(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(z1d);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= Xtol are removed. NaN
%  no data points in the gaps larger than Xtol are kept.

x1dv=X1D;x1dv(isnan(z1d) & sz <= Xtol)=[];
z1dv=z1d;z1dv(isnan(z1d) & sz <= Xtol)=[];

% interpolate to fill the small gaps and save as final subtransect profile

Z1Dsubtran(ntran,:)=interp1(x1dv,z1dv,X1D);

end % end subtransect loop

% calculate subtransect profile stats
Z1Dmean(nsurv,:)=mean(Z1Dsubtran,'omitnan');
Z1Dmed(nsurv,:)=median(Z1Dsubtran,'omitnan');
Z1Dstd(nsurv,:)=std(Z1Dsubtran,'omitnan');

end  % end survey loop

end




