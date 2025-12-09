% Build SM  Morpho struct files from SG grid files

mpath='/volumes/group/MOPS/';

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% input Mop number
%MopNumber=582;
%MopNumber=503;
for MopNumber=15:38%581:583%543:626%495:626
    %for MopNumber=582:582
    fprintf('%i of 626\n',MopNumber)

% load SG struct array
load([mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat' ],'SG');

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);

%-----------------------------------------------------------------
% ---- 2D parameter fields ----
%-----------------------------------------------------------------
% .X2D : x points of common gridded area for all surveys
% .Y2D : y points of common gridded area for all surveys
% .XY2Dnum : number of surveys with valid gridded data at a X2D,Y2D point 
% .Z2Dmin : x,y and minimum gridded z values using all surveys
% .Z2Dmax : x,y and maximum gridded z values using all surveys
% .Z2Dmean : mean gridded z surface, all surveys (mean of season means)

% get the min and max elevation grid points for all surveys

[ux, ~, xidx] = unique(vertcat(SG.X));
[uy, ~, yidx] = unique(vertcat(SG.Y));
z=vertcat(SG.Z);
% find largest class id (hardest substrate)
c=vertcat(SG.Class);
%c=vertcat(SG.X)*0;
cmax=accumarray([xidx(:), yidx(:)], c.',[],@nanmax);
% find min and max grid elevations
zcount = accumarray([xidx(:), yidx(:)], 1);  
zmin=accumarray([xidx(:), yidx(:)], z.',[],@nanmin);
zmax=accumarray([xidx(:), yidx(:)], z.',[],@nanmax);
ii=find(zcount > 0); % 1d indices of valid data
[i,j]=find(zcount > 0); % 2d indices of valid data
GM.Mopnum=MopNumber;
GM.X2D=ux(i);
GM.Y2D=uy(j);
GM.Z2Dmin=zmin(ii);
GM.Z2Dmax=zmax(ii);
GM.Z2Dclass=cmax(ii);
GM.XY2Dnum=zcount(ii);

% To calculate the mean elevation gridded points for all the surveys,
% first de-weight any experiment time periods of times of the year
% when there are larger numbers of surveys.  Month-of-year means are
% first calculated, and then the global mean is the mean on the month
% means.

for mnth=1:12
 GM.MM(mnth).X2D=[];
 GM.MM(mnth).Y2D=[];
 GM.MM(mnth).Z2Dmean=[];
 GM.MM(mnth).Z2Dmedian=[];
 midx=find(month(datetime(vertcat(SG.Datenum),'ConvertFrom','datenum')) == mnth);

 if ~isempty(midx)
 % bin and average rounded survey data by placing in unique
 %  x,y data array
  [ux, ~, xidx] = unique(vertcat(SG(midx).X));
  [uy, ~, yidx] = unique(vertcat(SG(midx).Y));
  z=vertcat(SG(midx).Z);
  %array of counts of the number of points at each unique x/y combination
  zcount = accumarray([xidx(:), yidx(:)], 1);  
  %array of average of z that fall into each unique x/y combination
  zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
  zmed=accumarray([xidx(:), yidx(:)], z.',[],@nanmedian);
  % reduce arrays to 1d vectors of x,y points with z data 
  ii=isnan(zavg(:)) == 0; % 1d indices of valid data
  [i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
  % final shore box data vectors
  % xutm=ux(i);yutm=uy(j);
  % zavg=zavg(ii);
  % add average data to monthly struct array
  GM.MM(mnth).X2D=ux(i);
  GM.MM(mnth).Y2D=uy(j);
  GM.MM(mnth).Z2Dmean=zavg(ii);
  GM.MM(mnth).Z2Dmedian=zmed(ii);
  
  % now calculate profile info for eeach monthly mean grid
  
  % reconstruct x,y grid points within gridded x,y area that 
  %  have "no data" NaN's

   xg=ux(i);yg=uy(j);zg=zavg(ii);
   nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
   tg=nan(ny,nx); % temp grid of Nans
   idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
   tg(idx)=zg; % add data to temp grid
   [X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);
   
   % now 2d interpolate z values for the Mop transect points
   zt=xt*NaN; % initialize transect elevation points 
   zt(:) = interp2(X,Y,tg,xt,yt); % interpolate transect points from temp grid
   GM.MM(mnth).X1D=x1d;
   GM.MM(mnth).Z1Dtransect=zt;
   
   % now 2d interpolate z values for all the subtransect points
   zst=xst*NaN; % initialize transect elevation points 
   zst(:) = interp2(X,Y,tg,xst(:),yst(:));
   % get mean, median, std, min-max z(xt) transect values for 20 subtransects.
   GM.MM(mnth).Z1Dmean=nanmean(zst,1);
   GM.MM(mnth).Z1Dmedian=nanmedian(zst,1);
   GM.MM(mnth).Z1Dstd=nanstd(zst,1);
   GM.MM(mnth).Z1Dmin=nanmin(zst,[],1);
   GM.MM(mnth).Z1Dmax=nanmax(zst,[],1);
 end
 
end

% now make mean surface based on monthly means
[ux, ~, xidx] = unique(vertcat(GM.MM.X2D));
[uy, ~, yidx] = unique(vertcat(GM.MM.Y2D));
z=vertcat(GM.MM.Z2Dmean);
zcount = accumarray([xidx(:), yidx(:)], 1);  
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
zstd=accumarray([xidx(:), yidx(:)], z.',[],@nanstd);
z=vertcat(GM.MM.Z2Dmedian);
zmed=accumarray([xidx(:), yidx(:)], z.',[],@nanmedian);
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% X2D=ux(i);
% Y2D=uy(j);
GM.Z2Dmean=zavg(ii);
GM.Z2Dmedian=zmed(ii);
GM.Z2Dstd=zstd(ii);
GM.XY2Dnummonths=zcount(ii);

% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

xg=GM.X2D;yg=GM.Y2D;zg=GM.Z2Dmean;cg=GM.Z2Dclass;
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);

% now 2d interpolate z values for the Mop transect points
zt=xt*NaN; % initialize transect elevation points 
zt(:) = interp2(X,Y,tg,xt,yt); % interpolate transect points from temp grid
  GM.X1D=x1d;
  GM.Z1Dtransect=zt;
% ct=xt*NaN; % initialize transect class points 
% ct(:) = interp2(X,Y,cg,xt,yt); % interpolate class points from temp grid
%   GM.Z1Dtransectclass=ct;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 20 subtransects.

  GM.Z1Dmean=nanmean(zst,1);
  GM.Z1Dmedian=nanmedian(zst,1);
  GM.Z1Dstd=nanstd(zst,1);
  GM.Z1Dmin=nanmin(zst,[],1);
  GM.Z1Dmax=nanmax(zst,[],1);

  tg(idx)=cg; % add class data to temp grid
% now 2d interpolate class values for all the subtransect points
  cst=xst*NaN; % initialize transect class points 
  cst(:) = interp2(X,Y,tg,xst(:),yst(:));
  % find largest class id (hardest substrate)
  GM.Z1Dclass=nanmax(cst,[],1);

% save global morpho results in GM file for Mop number
save(['M' num2str(MopNumber,'%5.5i') 'GM.mat'],'GM');

%-----------------------------------------------------------------
%  Individual survey morpho parmeters
%---------------------------------------------------------------

%-----------------------------------------------------------------
% ---- Volume fields ----
%-----------------------------------------------------------------
% .Xnear :  x points of common control area for nearshore volume estimates
% .Ynear :  y points of common control area for nearshore volume estimates
% .Vnear :  Volume (m^3) above Z2Dmin in common area 
% .Xbch :   x points defining area above MSL to the back beach line
% .Ybch :   y points defining area above MSL to the back beach line
% .Vbch :  Volume (m^3) above Z2Dmin and between MSL and the back beach


%-----------------------------------------------------------------
% ---- 1D parameter fields ----
%-----------------------------------------------------------------
% .X1D : xshore distance (m) from Mop back beach line
% .Z1Dtransect : gridded z interpolated on Mop transect from gridded data  
% .Z1Dmean : mean gridded z at xshore distance X1D
% .Z1Dmedian : mean gridded z at xshore distance X1D
% .Z1Dstd : standard deviation gridded z at xshore distance X1D
% .Z1Dmin : minimum gridded z at xshore distance X1D
% .Z1Dmax : minimum gridded z at xshore distance X1D

% loop through SG gridded data
clear SM

% load old SM struct array
load([mpath 'M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');
OM=SM;

for sn=1:size(SG,2)
    
    SM(sn).Mopnum=SG(sn).Mopnum;
    SM(sn).Datenum=SG(sn).Datenum;
    SM(sn).Source=SG(sn).Source;
    SM(sn).File=SG(sn).File;
    SM(sn).UTMzone=SG(sn).UTMzone;
    SM(sn).X1D=x1d;
    
    
    % Survey date and source already exist in old grid struct
    %  area, reassign to new SG struct array
    old=find( [OM.Datenum] == SM(sn).Datenum & ...
        strcmpi({OM.Source}, SM(sn).Source)==1);
  if isempty(old)  % new survey data, grid it and add to SG struct array
        fprintf(' New Data %s %s\n',datestr(SM(sn).Datenum),SM(sn).Source)
       
% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

if ~isempty(SG(sn).X)
xg=SG(sn).X;yg=SG(sn).Y;zg=SG(sn).Z;
cg=SG(sn).Class;
%cg=SG(sn).X*0;
%cmin=accumarray([xidx(:), yidx(:)], cg.',[],@nanmax);
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);

% now 2d interpolate z values for the Mop transect points
zt=xt*NaN; % initialize transect elevation points 
zt(:) = interp2(X,Y,tg,xt,yt);
  SM(sn).Z1Dtransect=zt;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 51 subtransects.

  SM(sn).Z1Dmean=nanmean(zst,1);
  SM(sn).Z1Dmedian=nanmedian(zst,1);
  SM(sn).Z1Dstd=nanstd(zst,1);
  SM(sn).Z1Dmin=nanmin(zst,[],1);
  SM(sn).Z1Dmax=nanmax(zst,[],1);
  
  tg(idx)=cg; % add class data to temp grid
  % now 2d interpolate class values for all the subtransect points
  cst=xst*NaN; % initialize transect class points 
  cst(:) = interp2(X,Y,tg,xst(:),yst(:));
  % find largest class id (hardest substrate)
  SM(sn).Z1Dclass=nanmax(cst,[],1);
end
  else
     SM(sn).Z1Dtransect=OM(old).Z1Dtransect;
     SM(sn).Z1Dmean=OM(old).Z1Dmean;
     SM(sn).Z1Dmedian=OM(old).Z1Dmedian;
     SM(sn).Z1Dstd=OM(old).Z1Dstd;
     SM(sn).Z1Dmin=OM(old).Z1Dmin;
     SM(sn).Z1Dmax=OM(old).Z1Dmax;
     SM(sn).Z1Dclass=OM(old).Z1Dclass;
  end
end
   
% save survey morpho results in SM file for Mop number
save([mpath 'M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM');

end

%--------------------------------------------------------------------
%--------------------------------------------------------------------
% 
% function [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,NumTrans,ExtendLine)
% 
% %  Calculates the interpolated 1m xshore resolution transect line, and NumTrans 
% %   subtransect lines, for/within the MopNumber survey area bounds. Extends
% %   the lines -/+ ExtendLine(1)/ ExtendLine(2) meters from the 
% %   back beach (-) and offshore (+) mop line, using the Mop table for
% %   reference.
% 
% % eg. for 20 subtransects that extend -100m back from back beach
% %      mop line:
% %
% %  [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,580,20,[-100 0])
% 
% %  returns transects all of the same length = longest of the
% %  subtransects
% 
% 
% %-------------------------------------------------------------
% % first, the subtransect lines
% %-------------------------------------------------------------
% 
%  % back beach neighboring mop midpoints
%  xb1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
%  yb1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
%  xb2=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
%  yb2=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
%  % offshore neighboring mop midpoints
%  xo1=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
%  yo1=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);
%  xo2=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
%  yo2=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
%  
%  nst=NumTrans; % number of subtransects 
%  % increase nst by 2 when calculating x,y step size to keep
%  %  first and last subtransects inside mop area boundaries
%  %  for 2d interpolation. 
%  xbst=linspace(xb1,xb2,nst+2);xbst=xbst(2:end-1);
%  ybst=linspace(yb1,yb2,nst+2);ybst=ybst(2:end-1);
%  xost=linspace(xo1,xo2,nst+2);xost=xost(2:end-1);
%  yost=linspace(yo1,yo2,nst+2);yost=yost(2:end-1);
%  
%  % get max subtransect length
%  dmax=0;
%  for n=1:nst
%      dist=pdist([xbst(n),ybst(n);xost(n),yost(n)]);
%      dmax=max([dmax dist]);
%  end
%  
%  % include the actual Mop transect distance 
%  x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
%  x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
%  dist=pdist([x1,y1;x2,y2]);
%  
%  dmax=max([dmax dist]); % max of the transect & subtransects lengths
%  dmax=ceil(dmax); % round up to whole meter
%  
% 
%  % now loop through subtransect endpoints and make 1m xshore interpolated
%  %  transects
%  for n=1:nst
%    x1=xbst(n);y1=ybst(n);x2=xost(n);y2=yost(n);
%    ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
%    dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
%    % adjust offshore extension so all transects have same xshore length
%    ExtendOff=ExtendLine(2)+(dmax-dist);
%  
%    % 1m spaced points along transect 
%    % x coords of 1m spaced points along transect
%    xst(n,:)=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
%    % y coords of 1m spaced points along transect
%    yst(n,:)=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
%  
%  end
%  
%  % xshore distance m from back beach point
%  x1d=ExtendLine(1):dmax;
%  
% %-------------------------------------------------------------
% % now the single Mop transect line
% %-------------------------------------------------------------
%     
%  % Mop transect 1m spaced xshore points
%  x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
%  x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
%  ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
%  dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
%  % adjust offshore extension so all transects have same xshore length
%  ExtendOff=ExtendLine(2)+(dmax-dist);
%  
%  % x coords of 1m spaced points along transect
%  xt=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
%  % y coords of 1m spaced points along transect
%  yt=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
%  
% end
