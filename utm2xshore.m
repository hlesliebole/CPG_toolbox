function [Xc,Zc]=utm2xshore(MopNumber,Xutm,Yutm)

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% divide mop area into 100 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,100,[-100 0]);

idx=find(vertcat(SA.Class) > 0);
Xutm=vertcat(SA.X);Yutm=vertcat(SA.Y);
Xutm=Xutm(idx);Yutm=Yutm(idx);
Z=vertcat(SA.Z);Z=Z(idx);

[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

[row,col] = ind2sub(size(xst),NearIdx);

hold on;plot(x1d(col),Z,'k.')
end