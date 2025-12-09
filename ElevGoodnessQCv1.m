function g=ElevGoodnessQCv1(x,y,z,dz)

% Version 1.  eg. g=ElevGoodnessQCv1(x,y,z,0.01)

% QC method for high spatial resolution LiDAR survey data.

% Assigns an elevation goodness score (g) to each x,y,z point based on where 
%  the z value falls in an elevation histogram compared to other points 
%  in the same survey elevation band polygon. Elevation band polygons are
%  defined by the x,y line bounding all the points that fall within a
%  a specified elevation range (but the resulting polygon can contain
%  points with elevations below/above this range). Contiguous elevation
%  ranges are used between the min and max survey elevations, with each pair 
%  of range limits chosen to insure that the resulting polygon will include 
%  at least N points (but may have many more, eg. in flat terrace areas).
%  The resulting polygons can overlap and the minimum goodness scores are
%  saved for points that fall within multiple polygons.
%
%  In a final step, the points that define the boundary line of the polygon 
%  that contains all the survey points is screened for outliers based on 
%  the difference in elevation between consecutive line points. 

% Input: x(n),y(n),z(n) data vectors

%   dz = histogram elevation "assessment" bin size (m). The smaller the bin
%        size, the more stringent the qc standard and the lower the overall 
%        goodness scores will be. Something in the neighborhood of
%        0.01 (1 cm) worked well with test 1m spatially averaged
%        x,y,z data from truck LiDAR in the tif files.  

% Output:  g(n) Elevation goodness factor for each x(n),y(n),z(n) point
%
%          g=0 is the lowest goodness and means a polygon elevation
%          histogram bin with n = 0 was encountered between the peak bin 
%          and the bin that includes the x,y,z point. ie. the point
%          is "separated" from the majority of the points in the polygon
%          by an "elevation air space" >= dz.
%
%          g=1 means a histogram bin with n = 1 was encountered between 
%          the peak bin and the bin that includes the x,y,z point. etc.
%          High g scores mean the point elevation is closer to the most
%          common elevations in the polygon.

%-------------------------------------------------------------------
 
%          The general qc strategy is to look for z outliers 
%          in a polygon that are above/below the min-max z range used to 
%          define the polygon area boundary.  So, it is assumed you
%          want to use min-max z ranges that are an order of magnitude
%          larger than the dz bin size used in the histograms to define
%          what is a significant incremental z difference.

polyZrange=(max(z)-min(z))/(dz*10);

% Because the beach profile points are not equally distrubuted by
%  elevation, convert the target z range to a minimum allowable
%  number of points in a polygon, N.  This will control the actual
%  (variable) z range step sizes as the algorithm progresses through
%  the x,y,z data from the min z to the max z.

% number of points in each polygon if the survey points were evenly
%  distributed across all z values.
N=round(length(z)/polyZrange); 

% initialize elevation goodness score array
gs=nan(size(z));

% sort data into ascending z order
[zs,isrt]=sort(z);xs=x(isrt);ys=y(isrt);

% number of elevation range steps
npts=length(zs);nsteps=floor(npts/N);Ns=ceil(npts/nsteps); 

% step through the sorted elevation points based on Ns
for n=1:Ns:nsteps*Ns  
    % min-max elevation of step, adjust end of last step to fit data length
    zmin=zs(n);if(n+Ns <= npts);zmax=zs(n+Ns);else;zmax=zs(end);end 
    idx=find(zs >= zmin & zs <= zmax); % point indices for this step
    k=boundary(xs(idx),ys(idx)); % indices of points on boundary of z range
    in=inpolygon(xs,ys,xs(idx(k)),ys(idx(k))); % out/in polygon logical
    % histogram bin elevations of points in polygon with dz resolution
    ipg=find(in == 1);[np,~,binnum]=histcounts(zs(ipg),...
         dz*floor(min(zs(ipg))/dz):dz:dz*ceil(100*max(zs(ipg))/dz));    
    % define goodness as cumulative bin minimums moving away from peak 
    [~,npk]=max(np); % find the peak bin 
    minlow=fliplr(cummin(fliplr(np(1:npk))));minhigh=cummin(np(npk:end));
    bgood=[minlow minhigh(2:end)]; % elevation goodness 
    % assign goodness scores to the in-polygon points. Keep minimum
    %  score for points that end up in multiple polygons.
    gs(ipg)=nanmin(horzcat(gs(ipg),bgood(binnum)'),[],2);
end

% finally, screen the boundary line of the polygon that encloses 
%  the entire survey for outlier points based on their boundary slopes
%  rather than elevations.
k=boundary(xs,ys);
[np,~,binnum]=histcounts(diff([zs(k)' zs(1)]));
[~,npk]=max(np);
devlow=fliplr(cummin(fliplr(np(1:npk))));devhigh=cummin(np(npk:end));
bgood=[devlow devhigh(2:end)]; % slope goodness 
gs(k)=nanmin(horzcat(gs(k),bgood(binnum)'),[],2);

g(isrt)=gs; % reassign goodness scores to the input unsorted indices 

end
