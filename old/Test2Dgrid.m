
% convert gridded x,y,z points to 2d matrices
x=GM.X2D;
y=GM.Y2D;
z=GM.Z2Dmin;

[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;
Z(idx)=z;
figure
surf(X,Y,Z,'linestyle','none','facecolor',[.5 .5 .5],'facealpha',.9);
hold on;
Z(idx)=GM.Z2Dmean;
surf(X,Y,Z,'linestyle','none','facecolor','b','facealpha',.5);
Z(idx)=GM.Z2Dmedian;
surf(X,Y,Z,'linestyle','none','facecolor','g','facealpha',.5);
Z(idx)=GM.Z2Dmax;
surf(X,Y,Z,'linestyle','none','facecolor','r','facealpha',.5);
legend('Min','Mean','Median','Max','location','northwest')

