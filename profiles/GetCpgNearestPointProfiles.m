function [X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol)

% usage:
%
%[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
%  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol)

% Dependencies:
%
%  needs /MOPS and MOPS/toolbox in its path
%  calls GetTransectLines.m and gapsize.m
%  uses  MopTableUTM.mat


%% --------------------------------------------------------------------

%  For each survey in the SA struct array, calculates mean and median 
%  profiles, and xshore standard deviations of a series of subtransect 
%  profiles spanning a Mop's area. 
%
%  The Mop area is divided into NumSubTrans subtransects centered on the
%  main Mop transect line. For each survey, nearest points-based profiles 
%  are calculated for each of the subtransects, and the mean,
%  median and standard deviation of the Mop area's subtransect xshore 
%  elevations are derived.
%

%% --------------------------------------------------------------------

% Input:

%   SA = 1m spatial avg struct array data file containing size(SA,2) surveys.
%   NumSubTrans = number of suntransects to divide the Mop area into
%   XgapTol = max subtransect xshore dist tolerance for small gap filling interpolation (m)
%   YdistTol = alongshore distance from mop subtransect tolerance (m)


%  Returns:
% X1Dmop(m) = Profile x indices (1m spaced integer values, x=0 is Mop Back Beach Point)
% X1Dcpg(m) = Profile x indices (1m spaced integer values, x=0 is CPG LiDAR Beach Only 
%                 data set boundary)

% Zdatetime(n) = MATLAB datetime variable of n survey dates
% Z1Dtrans(n,m) = MOP along-transect cross-shore elevations (meters, NAVD88)
% Z1Dmean(n,m) = MOP area mean cross-shore elevations
% Z1Dmedian(n,m) = MOP area median cross-shore elevations
% Z1Dmin(n,m) = MOP area cross-shore elevation minimums
% Z1Dmax(n,m) = MOP area cross-shore elevation maximums
% Z1Dstd(n,m) = MOP area cross-shore elevation standard deviations

%% --------------------------------------------------------------------

%% 1. Define the Mop Area subtransects using GetTransectLines. 
%   The subtransects follow the trapzoidal shape of Mop areas in an
%   accordian-like way.

% make NumSubTrans the next largest odd number if an even number was input
%  to the function. This makes the middle subtransect fall on the Mop line.
NumSubTrans=2*(floor(NumSubTrans/2))+1;

%  load Mop Transect Info in matlab tab;e varoable "Mop"
load('MopTableUTM.mat','Mop');
 
%   Get the main transect as UTM xt,yt points and subtransect lines as 
%   UTM xst,yst points ~1m apart.  This also sets the fixed 1m spaced X1D transect 
%   points starting -100m behind the mop back beach point out to the offshore point
[X1D,xt,yt,xst,yst]=GetTransectLines(Mop,SA(1).Mopnum,NumSubTrans,[-100 0]);

%% ----------------------------------------------------------
%% 2. Calculate nearest point profiles on subtransect lines

% initialize 2d output variables
Z1Dmean=NaN(size(SA,2),numel(X1D));
Z1Dtrans=Z1Dmean;
Z1Dmedian=Z1Dmean;
Z1Dstd=Z1Dmean;
Z1Dmin=Z1Dmean;
Z1Dmax=Z1Dmean;


% loop through all surveys
fprintf('\nCalculating profiles for %i Surveys using %i subtransects on each...\n',size(SA,2),NumSubTrans)

for nsurv=1:size(SA,2)

    % reinitialize the subtransect profile matrix for each survey 
    Z1Dsubtran=NaN(NumSubTrans,numel(X1D));

    %  Calculate a nearest point profile for each subtransect
    %
    % The nearest points profile variables are
    %
    %     Xuniq = the unique xshore distances with survey data (1m integers)
    %     Znear = the elevation of the survey point nearest to each Xuniq point 
    %     Zdist = the distance of each Znear point from the mop line (m)
    %
    %   Generates mop line profile vectors x1d and z1di that use a 
    %   distance from the mop line threshhold, YdistTol, to define valid profile 
    %   Znear data, and a maximum allowable spacing between Xuniq points, XgapTol, 
    %   that can be filled using 1d interpolation of the remaining Znear points,

    %% subtransect profile loop 
    for ntran=1:NumSubTrans
    clear Zdist Znear
    Z1Dsubtran(ntran,:)=NaN;
    
    % find the distances of all data apoints to the subtransect line points  
     [dp,NearIdx]=pdist2([yst(ntran,:)',xst(ntran,:)'],...
         [double(SA(nsurv).Y),double(SA(nsurv).X)],...
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
    % mop line (YdistTol) and xshore spatial gap (XgapTol) tolerances
    
    % set nearest points with Zdist > YdistTol to NaNs
    if exist('Zdist','var')
     Znear(Zdist > YdistTol)=NaN;
    
    
    % seed the regular spaced profile vectors (x1d and z1d) with the 
    %  valid nearest point elevations (both valid and NaNs)
     ndx=find(ismember(X1D,Xuniq));
     z1d=X1D*NaN;z1d(ndx)=Znear;
    
    % identify the size of the gap each xshore point is in (0 = not in a gap)
     sz=gapsize(z1d);
    
    % make vectors of valid data points to be used in interpolation
    %  where the x,z profile points in gaps <= XgapTol are removed. NaN
    %  no data points in the gaps larger than XgapTol are kept.
    
     x1dv=X1D;x1dv(isnan(z1d) & sz <= XgapTol)=[];
     z1dv=z1d;z1dv(isnan(z1d) & sz <= XgapTol)=[];
    
    % interpolate to fill the small gaps and save as final subtransect profile
    
     Z1Dsubtran(ntran,:)=interp1(x1dv,z1dv,X1D);
    
    % find the data back beach boundary for this subtransect
    xb=X1D(find(~isnan(Z1Dsubtran(ntran,:)), 1 ));
    if ~isempty(xb)
     xback(ntran)=xb;
    else
     xback(ntran)=NaN;
    end
    
    end
    
    end % end subtransect loop
    %% 


  % assign middle subtransect to the Mop transect line 
  mdl=median(1:NumSubTrans);
  Z1Dtrans(nsurv,:)=Z1Dsubtran(mdl,:);

  % define back beach as median location of substransects
  X1Dback(nsurv)=median(xback,'omitnan');

  % calculate subtransect profile stats if NumSubTrans
  if NumSubTrans > 1
      Z1Dmean(nsurv,:)=mean(Z1Dsubtran,'omitnan');
      Z1Dmedian(nsurv,:)=median(Z1Dsubtran,'omitnan');
      Z1Dstd(nsurv,:)=std(Z1Dsubtran,'omitnan');
      Z1Dmin(nsurv,:)=min(Z1Dsubtran);
      Z1Dmax(nsurv,:)=max(Z1Dsubtran);
  % otherwise use the transect result
  else
      Z1Dmean(nsurv,:)=Z1Dtrans(nsurv,:);
      Z1Dmedian(nsurv,:)=Z1Dtrans(nsurv,:);
      Z1Dstd(nsurv,:)=Z1Dtrans(nsurv,:)*0;
      Z1Dmin(nsurv,:)=Z1Dtrans(nsurv,:);
      Z1Dmax(nsurv,:)=Z1Dtrans(nsurv,:);
  end

end  % end survey loop

% use the median back beach data boundary from the trk and atv
% lidar surveys to define to lidar back beach and X1Dcpg
mobilelidar=find(strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR'));
 if ~isempty(mobilelidar)
  X0BeachOnly=median(X1Dback(mobilelidar),'omitnan');
 else
  X0BeachOnly=0;
 end

X1Dmop=X1D;
X1Dcpg=X1D-X0BeachOnly;

% return survey dates as matlab datetimes
Zdatetime=datetime([SA.Datenum],'ConvertFrom','datenum');

end




