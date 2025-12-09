function [Xm,Ym]=UTM2MopCoords(Xutm,Yutm)

Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);

% divide mop area into 101 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each. x1d is the xshore distance
%  from the backbeach line (line connecting Mop backbeach points).
%  xt,yt are the 2d grid locations of the main Mop transect line, and
%  xst,yst are the 2d grid locations for all the subtransects (including
%  the main line)
load MopTableUTM

Xm=Xutm*NaN;Ym=Yutm*NaN; % initialize mop coord vectors

for mop=unique(Nmop) % loop though range of mop numbers
    
 n=find(Nmop == mop); % points closest to this mop number

% get mpo subtransect lines with order 1m alongshore resolution
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,mop,101,[-100 500]);

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
  pdist2([yst(:),xst(:)],[double(Yutm(n)),double(Xutm(n))],'euclidean','smallest',1);

% define X based on the nearest transect line point
[row,col] = ind2sub(size(xst),NearIdx);
Xm(n)=x1d(col); % m cross-shore x
Ym(n)=mop-.5+(row-1)*0.01; % fractional Mop Y

end


end