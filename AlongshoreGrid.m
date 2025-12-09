function [Xa,Ya,Za]=AlongshoreGrid(Xutm,Yutm,Z)

% Alongshore gridding code
%
% 1. transform Mop area x,y utm (1m spatial res) grid points to
%    XshoreMop,Fractional Mop Number coordinates (1m xshore, 0.01 
%    Mop number resolution)
%
% 2. transform SA x,y utm points into XshoreMop,Fractional Mop Number
%    coordinates (also 1m xshore, 0.01 Mop NUmber res)
%
% 3. interpolate on the Mop number (alongshore) axis only to get area
%     grid point values.
%
% 4. save area grid points with valid data in their utm coords.
%
%
% redundant transformed x,y points are possible with compression
%  of alongshore coords on a curved coastline

% MopNumber=582;
% load M00582SA.mat
% ShoalsIdx=find(strcmp({SA.Source},'USACE'));
% JumboIdx=find(contains({SA.File},'umbo'));
% 
% n=JumboIdx(end);
% ndate=SA(JumboIdx(end)).Datenum;
% Xutm=SA(n).X;
% Yutm=SA(n).Y;
% Z=SA(n).Z;
% 
% % find same survey and adjacent mops
% load M00580SA.mat
% JumboIdx=find(contains({SA.File},'umbo'));
% idx=find([SA(JumboIdx).Datenum] == ndate);
% n=JumboIdx(idx);
% Xutm=vertcat(Xutm,SA(n).X);
% Yutm=vertcat(Yutm,SA(n).Y);
% Z=vertcat(Z,SA(n).Z);
% 
% % find same survey and adjacent mops
% load M00581SA.mat
% JumboIdx=find(contains({SA.File},'umbo'));
% idx=find([SA(JumboIdx).Datenum] == ndate);
% n=JumboIdx(idx);
% Xutm=vertcat(Xutm,SA(n).X);
% Yutm=vertcat(Yutm,SA(n).Y);
% Z=vertcat(Z,SA(n).Z);
% 
% % find same survey and adjacent mops
% load M00583SA.mat
% JumboIdx=find(contains({SA.File},'umbo'));
% idx=find([SA(JumboIdx).Datenum] == ndate);
% n=JumboIdx(idx);
% Xutm=vertcat(Xutm,SA(n).X);
% Yutm=vertcat(Yutm,SA(n).Y);
% Z=vertcat(Z,SA(n).Z);
% 
% load M00584SA.mat
% JumboIdx=find(contains({SA.File},'umbo'));
% idx=find([SA(JumboIdx).Datenum] == ndate);
% n=JumboIdx(idx);
% Xutm=vertcat(Xutm,SA(n).X);
% Yutm=vertcat(Yutm,SA(n).Y);
% Z=vertcat(Z,SA(n).Z);

%load MopTableUTM

[Xm,Ym]=UTM2MopCoords(Xutm,Yutm);

% check for duplicate x, y, z
     [uxy,ia,ic]=unique([Xm,Ym],'rows');uniqxyz=size(uxy,1);
     % remove those points
     Xm=Xm(ia);Ym=Ym(ia);Z=Z(ia);

%figure;plot(Xm,Ym,'k+')

XshoreWin=10;
AshoreWin=1;
xmin=floor(min(Xm));xmax=ceil(max(Xm));
ymin=floor(min(Ym*100));ymax=ceil(max(Ym*100));
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);
Zg=X*nan;

yq=0.01*(ymin:ymax);
n=0;
for x=xmin:xmax
    n=n+1;
    xq=x*ones(1,numel(yq));
    idx=find(Xm >= x-XshoreWin & Xm <= x+XshoreWin);
    % skip if collinear
    if sum(diff(Xm(idx))) > 0 
      zq=griddata(Xm(idx),Ym(idx),Z(idx),xq,yq);
      if numel(zq) == numel(xq)
          Zg(:,n)=zq;
      end
    end
end

%figure;imagesc(Zg);demcmap(Zg);colorbar;set(gca,'ydir','normal')  

Zg2=Zg;
xq=xmin:xmax;
n=0;
for y=ymin:ymax
    n=n+1;
    yq=y*ones(1,numel(xq));
    idx=find(Y(:) >= y-AshoreWin & Y(:) <= y+AshoreWin & ~isnan(Zg(:)));
    % skip if only collinear points
    if sum(diff(Y(idx))) > 0 
        zq=griddata(X(idx),Y(idx),Zg(idx),xq,yq);
        if numel(zq) == numel(xq)
          Zg2(n,:)=zq;
        end
    end
end

%figure;imagesc(Zg2);demcmap(Zg2);colorbar;set(gca,'ydir','normal')    
    
% convert gridded data back to utm coords

idx=find(~isnan(Zg2(:)));
[Xa,Ya]=Mop2UTMcoords(X(idx),Y(idx)/100);
Za=Zg2(idx);

end

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

function [Xutm,Yutm]=Mop2UTMcoords(Xm,Ym)

 Nmop=round(Ym); % Mops with data points

% divide mop area into 101 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each. x1d is the xshore distance
%  from the backbeach line (line connecting Mop backbeach points).
%  xt,yt are the 2d grid locations of the main Mop transect line, and
%  xst,yst are the 2d grid locations for all the subtransects (including
%  the main line)
 load MopTableUTM

 Xutm=Xm*NaN;Yutm=Ym*NaN; % initialize utm coord vectors

for mop=unique(Nmop)' % loop though range of mop numbers

% get mop subtransect lines with order 1m alongshore resolution
 [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,mop,101,[-100 500]);
    
 n=find(Nmop == mop); % index of points closest to this mop number

% x,y mop coords to subtransect line xst,yst utm location
 ix1d=round(Xm(n))-x1d(1)+1; % Mop coords xshore indices
 ix1d(ix1d < 1)=1;ix1d(ix1d > size(xst,2))=size(xst,2); % bounds check
 ist=round((Ym(n)-mop+.5)*100)+1; % subtransect indices 
 ist(ist < 1)=1;ist(ist > size(xst,1))=size(xst,1); % bounds check
 idx=sub2ind(size(xst),ist,ix1d); % 1d indices of subtransect points
% utm coords 
 Xutm(n)=xst(idx); 
 Yutm(n)=yst(idx);

end

end   




