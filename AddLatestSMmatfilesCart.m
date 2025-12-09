% Build SM  Morpho struct files from SG grid files

% folder containing CPG Mop mat files

MopPath='/volumes/group/MOPS/';

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% input Mop number
%MopNumber=582;
%MopNumber=503;
for MopNumber=578:584%495:626
    
% load SM struct array
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');
% load SG struct array
load(['M' num2str(MopNumber,'%5.5i') 'SG.mat' ],'SG');

newsrv=size(SG,2)-size(SM,2);
if newsrv > 0

fprintf('%i : adding %i new gridded surveys\n',MopNumber,newsrv)
    
% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);


%-----------------------------------------------------------------
%  Individual survey morpho parameters
%---------------------------------------------------------------


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

% loop through new SG gridded data

for sn=size(SM,2)+1:size(SG,2)
    fprintf('Survey %i\n',sn)
    SM(sn).Mopnum=SG(sn).Mopnum;
    SM(sn).Datenum=SG(sn).Datenum;
    SM(sn).Source=SG(sn).Source;
    SM(sn).File=SG(sn).File;
    SM(sn).UTMzone=SG(sn).UTMzone;
    SM(sn).X1D=x1d;
    
% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

if ~isempty(SG(sn).X)
xg=SG(sn).X;yg=SG(sn).Y;zg=SG(sn).Z;cg=SG(sn).Class;
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
else
    fprintf('SG survey X data was empty.\n',sn)
end
   
end
   
% save survey morpho results in SM file for Mop number
save([MopPath 'M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM');

else
    fprintf('%i : No new gridded surveys\n',MopNumber)
end
end

