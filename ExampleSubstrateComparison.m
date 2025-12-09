% some example code to compare substrate data sets
%  using the shoals 2009 and 2014 surveys as separate
%  estimates rather than as the combo dataset.

% load Shoals 2009 and 2014 roughness values
%  for just Mop 700
%  X,Y are utm eastings and northing 2d arrays
%  stdZ1 is 2009 roughness
%  stdZ2 is 2014 roughness
MopNumber=700;
ExampleRoughnessShoalsComparison

% cutoff between smooth (sand-gravel & low relief bedrock)
%  and rocky substrate
rcut=0.25;

% identify and combine rocky shoals points from 2009 and 2014
%  surveys
load Shoals2009and2014roughnessGrids.mat
% toss smooth values from both surveys
stdZ1(stdZ1 < rcut)=0;stdZ1(isnan(stdZ1))=0;
stdZ2(stdZ2 < rcut)=0;stdZ2(isnan(stdZ2))=0;
% make combined data set
stdZ3=stdZ1*0;
% at any shared 1m x,y resolution points keep
%  the max observed roughness above rcut
stdZ3(:)=max(vertcat(stdZ1(:)',stdZ2(:)'));
stdZ3(stdZ3 == 0)=NaN;
figure;pcolor(stdZ3);shading flat;colormap(jet)

% grid rocky points with max gap space filling of 10m
idx=find(~isnan(stdZ3(:)));
MaxGap=10; % set max gap between points for gridding at 250m
 [x,y,z]=addNoDataAreaPoints(X(idx),Y(idx),stdZ3(idx),MaxGap);         
% Grid the elevation shoals using Delaunay tesselation 
 zg=griddata(double(x),double(y),double(z),...
        double(min(X(:)):max(X(:))),double(min(Y(:)):max(Y(:)))');
% zg is gridded rocky areas with everything else NaNs
figure;surf(X,Y,zg);shading flat;colormap(jet);colorbar;view(2)

% now go back add smooth points around gridded rocky areas 
load Shoals2009and2014roughnessGrids.mat
% this time toss rocky points
stdZ1(stdZ1 >= rcut)=0;stdZ1(isnan(stdZ1))=0;
stdZ2(stdZ2 >= rcut)=0;stdZ2(isnan(stdZ2))=0;
stdZ4=stdZ1*0;
% at any shared 1m x,y resolution points keep
%  the max observed smooth value below rcut
stdZ4(:)=max(vertcat(stdZ1(:)',stdZ2(:)'));
stdZ4(stdZ4 == 0)=NaN;
figure;pcolor(stdZ4);shading flat;colormap(jet)
stdZ4(isnan(stdZ4))=0;

% now blend these smooth points with the gridded rocky data
zg(isnan(zg))=0;
% at any shared 1m x,y resolution points keep
%  the max value (the rocky value) 
stdZ4(:)=max(vertcat(stdZ4(:)',zg(:)'));
stdZ4(stdZ4 == 0)=NaN;
% the result is a combo of gridded rocky info with the 
%  no data gaps now containing smooth survey points
figure;pcolor(stdZ4);shading flat;colormap(jet);colorbar

% finally grid this combo data set with the same 10m gap filling
%  limit for the smooth points
idx=find(~isnan(stdZ4(:)));
MaxGap=10; % set max gap between points for gridding at 10m
 [x,y,z]=addNoDataAreaPoints(X(idx),Y(idx),stdZ4(idx),MaxGap);         
% Grid the elevation shoals using Delaunay tesselation 
 zg2=griddata(double(x),double(y),double(z),...
        double(min(X(:)):max(X(:))),double(min(Y(:)):max(Y(:)))');

% the result is gridded smooth (< rcut) and rocky (> rcut) values 
%   areas > 10m in size where there were no survey points in either
%   the 2009 & 2014 surveys have NaNs.  Some of these are likely
%   kelp canopy or dense algae (depths > few m's) or surfzone dropouts
%   (depth < few meters).
figure;surf(X,Y,zg2);shading flat;colormap(jet);colorbar;view(2)



% stdZ3=stdZ1+stdZ2;stdZ3(stdZ3 == 0)=NaN;
% 
% % grid stdZ1 data, use limits of X and Y for shared grid dimensions
% %   with stdZ2
% 
% idx=find(~isnan(stdZ1(:)));
% % bound the large data gaps with NaNs to avoid gridding them at this stage
%  MaxGap=250; % set max gap between points for gridding at 250m
%  [x,y,z]=addNoDataAreaPoints(X(idx),Y(idx),stdZ1(idx),MaxGap);         
% % Grid the elevation shoals using Delaunay tesselation 
%  zg1=griddata(double(x),double(y),double(z),...
%         double(min(X(:)):max(X(:))),double(min(Y(:)):max(Y(:)))');
%  
%  % grid stdZ2 data using same X Y grid dimensions
%  idx=find(~isnan(stdZ2(:)));   
%  [x,y,z]=addNoDataAreaPoints(X(idx),Y(idx),stdZ2(idx),MaxGap);   
%  zg2=griddata(double(x),double(y),double(z),...
%         double(min(X(:)):max(X(:))),double(min(Y(:)):max(Y(:)))');
%     
%  figure;surf(X,Y,zg1);shading flat;colormap(jet);colorbar;
%  set(gca,'clim',[0 .5],'xlim',[min(X(:)) max(X(:))],...
%      'ylim',[min(Y(:)) max(Y(:))]);view(2)
%  title(['Mop ' num2str(MopNumber) '  - Shoals 2009'])
%  figure;surf(X,Y,zg2);shading flat;colormap(jet);colorbar;
%  set(gca,'clim',[0 .5],'xlim',[min(X(:)) max(X(:))],...
%      'ylim',[min(Y(:)) max(Y(:))]);view(2)
%  title(['Mop ' num2str(MopNumber) '  - Shoals 2014'])
%    figure;surf(X,Y,zg2-zg1);shading flat;colormap(jet);colorbar;
%  set(gca,'xlim',[min(X(:)) max(X(:))],...
%      'ylim',[min(Y(:)) max(Y(:))]);view(2)
%  title(['Mop ' num2str(MopNumber) '  - Shoals 2014 minus 2009']) 
%  
%  figure;histogram2(zg1(:),zg2(:),'displaystyle','tile')
%  
%  idx=find(zg1(:) > 0.2 & zg2(:) > 0.2);
%  zr=zg1*NaN;zr(idx)=max(vertcat(zg2(idx)',zg1(idx)'));
%  figure;surf(X,Y,zr);shading flat;colormap(jet);colorbar;
%   set(gca,'clim',[0 .5],'xlim',[min(X(:)) max(X(:))],...
%      'ylim',[min(Y(:)) max(Y(:))]);view(2)
%  
% %      [YG,XG]=find(~isnan(zg)); % find valid grid data x,y vectors
% %      Vzg=zg(~isnan(zg(:))); % valid zg data vector
% %      XGutm=min(x)-1+XG; % adjust back to utm coords
% %      YGutm=min(y)-1+YG;
% 
