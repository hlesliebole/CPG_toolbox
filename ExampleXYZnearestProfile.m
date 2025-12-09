% Example code to "project" UTM x,y,z data onto a transect defined by
% a back beach UTM point xb,yb and an offshore UTM point xoff,yoff.
% 
% The result is 2 vectors: X1D and Z1D defining the transect profile
%  with 1m xshore resolution.

% The method recipe is:
%   1. define a discrete transect of UTM xt,yt points that are 1m apart
%        along the transect line, X1D(-100:1:N).  X1D(0) represents xb,yb,
%        and X1D(N) represents xoff,yoff , N (rounded to integer) meters away
%        from xb,yb
%   2. "Project" nearest x,y,z elevation points, no more that Ytol away, into 
%        the discrete 1m xshore bins of X1D and Z1D.  If multiple points fall 
%        in the same xshore bin, the nearest one is used.
%   3. Fill any X1D data gaps less than or equal to Xtol using linear
%        interpolation of Z1D.

%% --- Settings ---

Ytol=3; % Max distance (m) a point can be from the transect line
Xtol=3; % Max xshore gap size (m) to be filled by interpolation. Otherwise
        %  the Z1D value will be a NaN.

%% define transect end points in UTM
%  for this example use Mop 511 by the SIO lifeguard tower
load('MopTableUTM.mat','Mop');
xb=Mop.BackXutm(511);yb=Mop.BackYutm(511);
xoff=Mop.OffXutm(511);yoff=Mop.OffYutm(511);

%% define x,y,z data points in UTM
%  use the last truck lidar survey for Mop 511 
%  saved in group/MOPS/M00511SA.mat on reefbreak1
load('M00511SA.mat','SA');
Trkidx=find(contains({SA.Source}','Trk'));
x=SA(Trkidx(end)).X;
y=SA(Trkidx(end)).Y;
z=SA(Trkidx(end)).Z;

% ---------------
%%  1. Define the UTM xt,yt 1m transect points, and the X1D xshore distance 
%%     vector, based on input xb,yb and xoff,yoff
% ---------------

% option to extend the transect line landward and seaward of the input end
% points.  For Mops, I typically allow for up to 100m behind the back beach
% point but don't extend seaward.
ExtendLine=[-100 +0]; % use integer values

 % Define X1D transect xshore dimension
 Dist=ceil(pdist([xb,yb;xoff,yoff])); %round distance between end points to 1m
 X1D=ExtendLine(1):Dist+ExtendLine(2); % add extensions
 % define integer offshore end point  
 ExtendOff=ExtendLine(2)+Dist; 

 % Define transect points with 1m spacing
 ang=atan2(yoff-yb,xoff-xb); % radian angle between back and offshore point
 
 % xt UTM coords of 1m spaced points along transect
 xt=xb+ExtendLine(1)*cos(ang):cos(ang):xoff+ExtendOff*cos(ang);
 % yt UTM coords of 1m spaced points along transect
 yt=yb+ExtendLine(1)*sin(ang):sin(ang):yoff+ExtendOff*sin(ang);  

% ---------------
%% 2. Project data points onto the transect points
% ---------------

% rather than formally projecting points using the transect perpendicular, find the
%  closest transect point to each data point instead, which will represent 
%   the data points perpendicular projection

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yt(:),xt(:)],[y,x],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xt),NearIdx); 

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
    Znear(n)=z(idx(imin));   
end

% make an interpolated profile based on distance from 
% transectline (Ytol) and xshore spatial gap (Xtol) tolerances

% set nearest points with Zdist > dy to NaNs
Znear(Zdist > Ytol)=NaN;

% seed the regular spaced profile vectors (X1D and Z1D) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(X1D,Xuniq));
Z1D=X1D*NaN;Z1D(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(Z1D);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= Xtol are temporarily removed. 
%  NaN no data points in the gaps larger than Xtol are kept.

X1Dv=X1D;X1Dv(isnan(Z1D) & sz <= Xtol)=[];
Z1Dv=Z1D;Z1Dv(isnan(Z1D) & sz <= Xtol)=[];

% interpolate to fill the small gaps
Z1D=interp1(X1Dv,Z1Dv,X1D);

figure;plot(X1D,Z1D,'k*-');
title([SA(Trkidx(end)).Source ' ' datestr(SA(Trkidx(end)).Datenum)])

%% ---------------------------------------------------------------

function sz=gapsize(x)
% Calculates the gap size of a vector. 
%
% A gap is defined as the number of consequtive nans.
%
% USAGE: sz=gapsize(x)
% sz is same length as input. 
%
% example:
% x=rand(20,1);
% x(x>.5)=nan; 
% [x gapsize(x)]
%
% aslak grinsted 2010
x=~isnan(x);
hasdata=[0;find(x(:)); length(x)+1];
sz=zeros(size(x));
for ii=1:length(hasdata)-1
    ix=hasdata(ii)+1:hasdata(ii+1)-1;
    sz(ix)=length(ix);
end

end



