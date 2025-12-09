
% read in 10/7/2021 jetski survey
dfile='/volumes/group/ig8_wheel/20211011_00581_00583_Torrey_wheel/filtered_clean20211011.llnezts.navd88';

% read in data
[xutm,yutm,z,c]=readSurveyFileUTM('gps',dfile);

%----------------------------------------
% reduce to 1m spatial averages  
%----------------------------------------
 
Res=1; % 1m spatial resolution

% round survey x,y to desired resolution
xr=Res*round(xutm/Res); % round to Res meter spatial resolution
yr=Res*round(yutm/Res); % 

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(xr);
[uy, ~, yidx] = unique(yr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average of z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
%tavg = accumarray([xidx(:), yidx(:)], t.')./zcount;
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final shore box data vectors
xutm=ux(i);yutm=uy(j);
zavg=zavg(ii);
cmode=cmode(ii);
%tavg=tavg(ii);

%fprintf('Min-Max Classification Mode: %d %d\n',min(cmode),max(cmode))

%--------------------------------------------------------
% option to remove spatially averaged data based on less 
%  than N survey points
%
% zcount=zcount(ii); 
% N=3;
%i=find(zcount < N);xutm(i)=[];yutm(i)=[];zavg(i)=[];
%  zcount(i)=[];cmode(i)=[];
%--------------------------------------------------------

%end % end if for npts > 0

fprintf(1,...
    'Reduced to %g , %g x %g meter spatially averaged survey points\n',...
    length(xutm),Res,Res);

%-------------------------
% grid the aggregated data

Xr=xutm;Yr=yutm;Zr=zavg;

     fprintf('Gridding %i points...\n',length(Xr))
     
% bound the data gaps with NaNs to avoid gridding them
     MaxGap=250;
     fprintf('Max Data Gap to grid across = %d meters\n',MaxGap)
    
     [x,y,z]=addNoDataAreaPoints(Xr,Yr,Zr,MaxGap);     
    
% Grid the elevation survey using Delaunay tesselation 
    zg=griddata(double(x),double(y),double(z),...
        double(min(x):max(x)),double(min(y):max(y))');
    

figure
imagesc(min(x):max(x),min(y):max(y),flipud(zg))
colormap(jet);colorbar;axis equal