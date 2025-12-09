
% look at surface gradients of a gridded survey

%n=1;

allG=[];
allZ=[];
for n=1:size(SG,2)
    
xg=SG(n).X;yg=SG(n).Y;zg=SG(n).Z;

%xg=GM.X2D;yg=GM.Y2D;zg=GM.Z2Dmean;cg=GM.Z2Dclass;
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
Z=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
Z(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);
[fx,fy]=gradient(Z);
G=sqrt(fx.^2+fy.^2);
allG=[ allG G(:)' ];
allZ=[ allZ Z(:)' ];

end

histogram2(allG,allZ,'facecolor','flat')
set(gca,'clim',[0 50])
colormap(jet)
colorbar