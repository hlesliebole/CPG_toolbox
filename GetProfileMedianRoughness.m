function [x1d,Z1Dmedian,R1Dmedian]=GetProfileMedianRoughness(MopNumber)

% load SG struct array
%load([mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat' ],'SG');
%clearvars

%MopNumber=666;
WinSize=11;
MinWinPts=5;

RS=GetMopShoalsRoughnessV2(MopNumber,WinSize,MinWinPts);
sn=2;% 2014 only % combined shoals surveys roughness
if ~isempty(RS(sn).Xutm)
xg=RS(sn).Xutm;yg=RS(sn).Yutm;zg=RS(sn).Z;rg=RS(sn).Sigma;
% % reduce grid area to just include the data in mop range
xmin=min(xg);xmax=max(xg);
ymin=min(yg);ymax=max(yg);
[Xg,Yg]=meshgrid(xmin:xmax,ymin:ymax);
Zg=griddata(xg,yg,zg,xmin:xmax,[ymin:ymax]');%,'nearest');
Rg=griddata(xg,yg,rg,xmin:xmax,[ymin:ymax]');%,'nearest');

idx=find(~isnan(Zg(:)));
xg=Xg(idx);yg=Yg(idx);zg=Zg(idx);rg=Rg(idx);

load MopTableUTM.mat
% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);

nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
trg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
trg(idx)=rg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);


% now 2d interpolate z values for the Mop transect points
zt=xt*NaN; % initialize transect elevation points 
zt(:) = interp2(X,Y,tg,xt,yt);
  %SM(sn).Z1Dtransect=zt;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
rst=xst*NaN; % initialize transect roughness points 
rst(:) = interp2(X,Y,trg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 51 subtransects.

%  RS(sn).Z1Dmean=nanmean(zst,1);
  Z1Dmedian=nanmedian(zst,1);
  R1Dmedian=nanmedian(rst,1);
  % RS(sn).Z1Dstd=nanstd(zst,1);
  % RS(sn).Z1Dmin=nanmin(zst,[],1);
  % RS(sn).Z1Dmax=nanmax(zst,[],1);
  
end

end

