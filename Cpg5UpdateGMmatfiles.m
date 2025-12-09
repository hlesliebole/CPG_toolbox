% Build SM  Morpho struct files from SG grid files
clearvars
CpgDefineMopPath

load SurveyMasterListWithMops.mat % Load current Survey info
load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% find all the Mop numbers with survey data
MopNums=unique([Survey.NearestMops]);
MopNums=sort(MopNums);

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% input Mop number
%MopNumber=582;
%MopNumber=503;
for MopNumber=620:633%647:667%700:750%634:646%578:595%743:743%582:582%MopNums(29:end)%15:38%581:583%543:626%495:626
    
 fprintf('%i of %i\n',MopNumber,numel(MopNums))

if exist([mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat' ],'file')...
        && MopNumber > 1

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
 GM.MM(mnth).NumSurveys=0;
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
  GM.MM(mnth).NumSurveys=numel(midx);
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
save([mpath 'M' num2str(MopNumber,'%5.5i') 'GM.mat'],'GM');

end
end
