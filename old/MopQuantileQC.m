% function [DpZ,DpBeta,DpL]=MopBeachQuantileQC(x,y,z,dz)
dz=0.1;
load M00503SA.mat
n=6;
x=SA(n).X;y=SA(n).Y;z=SA(n).Z;
% Input: N, 1-m spatially averaged beach survey data points x(N), y(N), z(N)
%           in UTM coordinates, and the size of the discrete elevation bins
%           to use for QC.
%
% Output: quantile outlier probability difference DpZ(N),DpBeta(N),DpL(N)
%         for each point.

%----------------------------------------------------------
% Step 1: Put x,y,z spatially averaged data in a 2D elevation 
%         grid arrays, X,Y,Z. "no data" grid points = NaNs.
%----------------------------------------------------------

nx=max(x)-min(x)+3;ny=max(y)-min(y)+3; % grid dimensions
Z=nan(ny,nx); % initial grid of NaNs
gidx=sub2ind([ny nx],y-min(y)+2,x-min(x)+2); % data point 1d grid indices
Z(gidx)=z; % add elevation data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1); % 2D grid indices

%----------------------------------------------------------
% Step 2: Replace any NaNs within Z grid array's survey spatial 
%         bounds with interpolated elevations.
%----------------------------------------------------------

iData=find(~isnan(Z(:))); % 1d indices of grid points with data 
iNoData=find(isnan(Z(:))); % 1d indices of grid points without data
% delaunay interpolate to fill in any no data gaps in survey area
Z(iNoData)=griddata(X(iData),Y(iData),Z(iData),X(iNoData),Y(iNoData));

%-------------------------------------------------------------
% Step 3: Define extreme quantile exceedance differences DqZ
%         at all grid points.
%-------------------------------------------------------------

% get quantile function and min-max slope break probabilities
[qZ,p,pminz,pmaxz]=GetQuantCutoffs(Z,1); 
% interpolate quantile function p,qZ
pZ=nan(ny,nx);pZ(:)=interp1(qZ,p,Z(:));
% probability differences from min-max slope break points
DpZ=pZ;DpZ(:)=max([pminz-pZ(:)';pZ(:)'-pmaxz]);

%-------------------------------------------------------------
% Step 4: Temporarily replace any Z points that have DpZ > 0 with
%         NaN's and reinterpolate.
%-------------------------------------------------------------

iData=find(DpZ(:) <= 0); % 1d indices of grid points with good data 
iNoData=find(DpZ(:) > 0); % 1d indices of grid points with flagged data
% delaunay interpolate to fill in any no data gaps in survey area
Z(iNoData)=griddata(X(iData),Y(iData),Z(iData),X(iNoData),Y(iNoData));

%-------------------------------------------------------------
% Step 5: Loop through the survey's elevation bands and define
%         DpBeta at all grid points.
%-------------------------------------------------------------

yn=0;
for loop=1:5
    
[sx,sy]=gradient(Z);
%sx=movmax(abs(sx),21,1);sy=movmax(abs(sy),21,1);
sx=movmax(abs(sx),7,2);sy=movmax(abs(sy),7,1);
% use beta = max gradient component at each x,y location
B=max(cat(3,abs(sx),abs(sy)),[],3);
%B=movmax(B,loop,1);
%if(yn==1);B=movmax(sx,5,1);yn=1;else;B=movmax(sy,5,1);yn=0;end
% figure out elevation range of data rounded to nearest "dz" 
zmin=dz*floor(min(Z(:))/dz);zmax=dz*ceil(max(Z(:))/dz);
pB=nan(ny,nx);DpBeta=nan(ny,nx); % initialize slope array
for zl=zmin:dz:zmax-dz  % loop through elevation bins
    idx=find(Z(:) >= zl & Z(:) < zl+dz & ~isnan(B(:)));
    if ~isempty(idx) % skip if no elevs in this range
        [qB,p,pminB,pmaxB]=GetQuantCutoffs(B(idx),2);
        [qBu,ia,ic]=unique(qB);pu=p(ia); % remove duplicate function values
        %[Bu,ia,ic]=unique(B(idx));
        pB(idx)=interp1(qBu,pu,B(idx));
        %pB(idx)=pBsrt(ic);
        DpBeta(:)=max([pminB-pB(:)';pB(:)'-pmaxB]);
    end
end

length(find(DpBeta(:) > 0))

Z(DpBeta(:) > 0)=NaN;
iData=find(~isnan(Z(:))); % 1d indices of grid points with good data 
iNoData=find(isnan(Z(:))); % 1d indices of grid points with flagged data
% delaunay interpolate to fill in any no data gaps in survey area
Z(iNoData)=griddata(X(iData),Y(iData),Z(iData),X(iNoData),Y(iNoData),'linear');

end
% 
% %------------------------------------------------------------------
% % Step 6: Loop through the elevation bins between 1.0 and 1.5m NAVD88
% %         with at least 200 points and find the one with the best fit 
% %         (lowest error variance) to a straight line.  This line will 
% %         be used to define relative Xshore location distances for all 
% %         points.
% %-------------------------------------------------------------------
% 
% minstd=Inf;abest=[];bbest=[];zbest=[];
% zmin=dz*floor(1.0/dz);zmax=dz*ceil(1.5/dz);
% for zl=zmin:dz:zmax
%     idx=find(Z(:) >= zl & Z(:) < zl+dz); % points in elev range
%     if length(idx) > 200 % skip if too few elevs in this range
%         [a,b]=linfit(X(idx),Y(idx)); % best fit line to elv points
%         L=(a.*X(idx) - Y(idx) + b)./sqrt(a.^2+1); % distance from line
%         if minstd > nanstd(L) % standard dev sigma distance
%             minstd=nanstd(L);abest=a;bbest=b;zbest=zl+dz/2;     
%         end
%     end
% end
% 
% % best fit normal (True compass coming from) 
% BchNorm=360-atand(abest);
% 
% %--------------------------------------------------------------------
% % Step 7: Loop through the survey's elevation bands again and define
% %         DqL at all grid points.
% %---------------------------------------------------------------------
% 
% zmin=dz*floor(min(Z(:))/dz);zmax=dz*ceil(max(Z(:))/dz);
% pL=nan(ny,nx);DpL=nan(ny,nx); % initialize slope array
% for zl=zmin:dz:zmax-dz  % loop through elevation bins
%     idx=find(Z(:) >= zl & Z(:) < zl+dz & ~isnan(B(:)));
%     if ~isempty(idx) % skip if no elevs in this range
%       % distance from best line
%        L=(abest.*X(idx) - Y(idx) + bbest)./sqrt(abest.^2+1);
%        [qL,p,pminL,pmaxL]=GetQuantCutoffs(L,1);
%        [qLu,ia,ic]=unique(qL);pu=p(ia); % remove duplicate function values
%        pL(idx)=interp1(qLu,pu,L);
%        DpL(idx)=max([pminL-pL(idx)';pL(idx)'-pmaxL]);
%     end
% end
% 
% %--------------------------------------------------------------------
% % Step 8: Reduce the grids of extreme quantile exceedance differences 
% %         to vectors of values DpZ(N),DpBeta(N),DpL(N) at the N 
% %         spatially averaged input points.
% %---------------------------------------------------------------------
% 
% 
% 
% 
% 
% [nz,zn]=histcounts(Z(:),zmin:dz:zmax);
% zn=zn(1)+dz/2:dz:zn(end-1)+dz/2;
% 
% figure;
% plot(zn,nz,'k-',zn,nz,'k.','markersize',10);grid on
% xlabel('Bin Elevation (m, NAVD88)');
% ylabel('Number of Grid Points in Bin');
% hold on;
% 
% [dmin,dmax]=GetQuantCutoffs(Z,1);
% 
% plot([dmin dmin],[0 max(nz)],'k--');
% plot([dmax dmax],[0 max(nz)],'k--');
% 
% 
% [sx,sy]=gradient(Z);
% % use beta = max gradient component at each x,y location
% beta=max(cat(3,abs(sx),abs(sy)),[],3);
% 
% % loop through elevation levels and find the one
% %   with at least 200 points that has the best
% %   (smallest standard deviation) linear fit line.
% 
% minstd=Inf;abest=[];bbest=[];zbest=[];
% for z=zmin:dz:zmax-dz
%     idx=find(Z(:) >= z & Z(:) < z+dz); % points in elev range
%     if length(idx) > 200 % skip if too few elevs in this range
%         [a,b]=linfit(X(idx),Y(idx)); % best fit line to elv points
%         d=(a.*X(idx) - Y(idx) + b)./sqrt(a.^2+1); % distance from line
%         if minstd > nanstd(d) % standard dev sigma distance
%             minstd=nanstd(d);abest=a;bbest=b;zbest=z+dz/2;
%         end
%     end
% end
% 
% % best fit normal (True compass coming from) 
% BchNorm=360-atand(abest);
% 
% yyaxis right
% smed=[];
% for z=zmin:dz:zmax-dz
%     idx=find(Z(:) >= z & Z(:) < z+dz);
%     if ~isempty(idx) % skip if no elevs in this range
%     plot((z+dz/2)*ones(length(idx),1),s(idx),'g.');
%     plot((z+dz/2),nanmedian(s(idx)),'b.');
%     smed=[smed nanmedian(s(idx))];
%     [dmin,dmax]=GetQuantCutoffs(s(idx),2);
%     plot((z+dz/2),dmin,'r.');
%     plot((z+dz/2),dmax,'r.');
%     end
% end
% ylabel('Grid Point Gradients');
% %figure;plot(zn,smed)
