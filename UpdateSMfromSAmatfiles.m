% Build SM  Morpho struct files from SG grid files

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% input Mop number
%MopNumber=582;
%MopNumber=503;
for MopNumber=582:582%15:38%581:583%543:626%495:626
    %for MopNumber=582:582
    fprintf('%i of 626\n',MopNumber)

% load SG struct array
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat' ],'SA');

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);


%-----------------------------------------------------------------
%  Individual survey morpho parmeters
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

% loop through SG gridded data
clear SM

% load old SM struct array
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');
%OM=SM;

for sn=363:363%1:size(SA,2)
    
%     SM(sn).Mopnum=SG(sn).Mopnum;
%     SM(sn).Datenum=SG(sn).Datenum;
%     SM(sn).Source=SG(sn).Source;
%     SM(sn).File=SG(sn).File;
%     SM(sn).UTMzone=SG(sn).UTMzone;
%     SM(sn).X1D=x1d;
    
%     % Survey date and source already exist in old grid struct
%     %  area, reassign to new SG struct array
%     old=find( [OM.Datenum] == SM(sn).Datenum & ...
%         strcmpi({OM.Source}, SM(sn).Source)==1);
%   if isempty(old)  % new survey data, grid it and add to SG struct array
%         fprintf(' New Data %s %s\n',datestr(SM(sn).Datenum),SM(sn).Source)
       
% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

if ~isempty(SA(sn).X)
xg=SA(sn).X;yg=SA(sn).Y;zg=SA(sn).Z;
cg=SA(sn).Class;
%cg=SG(sn).X*0;
%cmin=accumarray([xidx(:), yidx(:)], cg.',[],@nanmax);
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);

% % now 2d interpolate z values for the Mop transect points
% zt=xt*NaN; % initialize transect elevation points 
% zt(:) = interp2(X,Y,tg,xt,yt);
%   SM(sn).Z1Dtransect=zt;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 51 subtransects.

  %SM(sn).Z1Drawmean=nanmean(zst,1);
  Z1Drawmean=nanmean(zst,1);
%   SM(sn).Z1Dmedian=nanmedian(zst,1);
%   SM(sn).Z1Dstd=nanstd(zst,1);
%   SM(sn).Z1Dmin=nanmin(zst,[],1);
%   SM(sn).Z1Dmax=nanmax(zst,[],1);
  
%   tg(idx)=cg; % add class data to temp grid
%   % now 2d interpolate class values for all the subtransect points
%   cst=xst*NaN; % initialize transect class points 
%   cst(:) = interp2(X,Y,tg,xst(:),yst(:));
%   % find largest class id (hardest substrate)
%   SM(sn).Z1Dclass=nanmax(cst,[],1);
end
%   else
%      SM(sn).Z1Dtransect=OM(old).Z1Dtransect;
%      SM(sn).Z1Dmean=OM(old).Z1Dmean;
%      SM(sn).Z1Dmedian=OM(old).Z1Dmedian;
%      SM(sn).Z1Dstd=OM(old).Z1Dstd;
%      SM(sn).Z1Dmin=OM(old).Z1Dmin;
%      SM(sn).Z1Dmax=OM(old).Z1Dmax;
%      SM(sn).Z1Dclass=OM(old).Z1Dclass;
%   end
end
   
% save survey morpho results in SM file for Mop number
%save(['M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM');

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
