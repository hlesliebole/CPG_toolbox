
% load M00503SA.mat
% n=1;
% xg=SA(n).X;yg=SA(n).Y;zg=SA(n).Z;
% look at surface gradients of a gridded survey

% allG=[];
% allZ=[];
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
Z=nan(ny,nx); % temp grid of Nans

idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
Z(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);
%[fx,fy]=gradient(Z);fx=fx(idx);fy=fy(idx);
dx1=diff(Z,1,2);dx1(:,end+1)=dx1(:,end); % x grid elev differences
dy1=diff(Z,1,1);dy1(end+1,:)=dy1(end,:); % y grid elevation differences
dx2=fliplr(diff(fliplr(Z),1,2));dx2(:,end+1)=dx2(:,end);
dy2=flipud(diff(flipud(Z),1,2));dy2(:,end+1)=dy2(:,end);
stpx=(abs(dx1)+abs(dx2))/2;
stpy=(abs(dy1)+abs(dy2))/2;
stp=sqrt(stpx.^2+stpy.^2);

Na=[];
S2=[];
S3=[];
S4=[];

for frac=.01:.01:1
%for N=[100:100:length(zg) length(zg)]
N=round(length(zg)*frac);    
%N=50; % target min sample size for each elevation step

StpDev=nan(ny,nx); % mean steepness deviation array
XshDev=nan(ny,nx); % xshore location deviation array



% round Z to nearest mm
Z=round(Z,3);
zmin=min(Z(:));zmax=max(Z(:));

z1=zmin;
z2=zmin;

ScatterPlotMopUTM(xg,yg,zg,'2d');
hold on;xl=get(gca,'xlim');

while z2 < zmax  % loop through elevation range
    z1=z2;
    Nmed=0;
while Nmed < N && z2 <= zmax  % find next elevation range that has >= N points
    z2=z2+0.001;
    idx=find(Z(:) >= z1 & Z(:) < z2);
    Nmed=length(idx);
end
   if(z2 > zmax && Nmed < N) % fix last sample bin
      zdiff=abs(Z(:)-max(Z(:))); % elev diff from point being assessed
      [zsort,zidx]=sort(zdiff); % sort from closest to furthest
      idx=zidx(1:N); % highest N elevation points
      Nmed=length(idx);
   %   fprintf('%5.3f %5.3f %5i\n',min(Z(idx)),max(Z(idx)),Nmed)
   else
   %  fprintf('%5.3f %5.3f %5i\n',z1,z2,Nmed)
   end
   
    %StpDev(idx)=stp(idx); 
   
    % Mean Steepness Deviation
    stpmean=nanmean(stp(idx)); % mean steepness of elev range
    stpsigma=nanstd(stp(idx)); % steepness standard dev sigma of elev range
    StpDev(idx)=(stp(idx)-stpmean)/stpsigma; % steepness deviations in fraction of sigma
      
    [a,b]=linfit(X(idx),Y(idx)); % best fit line to elv points
    p1=plot(xl,xl*a+b,'k-');hold on;
    d=(a.*X(idx) - Y(idx) + b)./sqrt(a.^2+1); % distance from line
    dmean=nanmean(d); % mean distance from best fit line
    dsigma=nanstd(d); % standard dev sigma distance
    XshDev(idx)=(d-dmean)/dsigma;
      
end
Na=[Na N];
S2=[S2 length(find(StpDev > 2))];
S3=[S3 length(find(StpDev > 3))];
S4=[S4 length(find(StpDev > 4))];


fprintf('%5i %10i %10i %10i\n',N,length(find(StpDev > 2)),...
   length(find(StpDev > 3)),length(find(StpDev > 4))) 
end

figure;pcolor(abs(StpDev));shading flat;colormap(jet);colorbar
set(gca,'clim',[0 3.5]);colormap(jet(7));

figure;pcolor(abs(XshDev));shading flat;colormap(jet);colorbar
set(gca,'clim',[0 3.5]);colormap(jet(7));

% poor mans fix to truncated normal distribution
figure; histogram([stp(:)' (2*nanmean(stp(:))-stp(ineg))'])

PERCENTRANK = @(YourArray, TheProbes) reshape( mean( bsxfun(@le, YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

% % figure out smallest step size that has a median sample
% %  size > N
% 
% dz=0;
% 
% zmin=min(Z(:));zmax=max(Z(:));
% 
% 
% while Nmed < N
%     dz=dz+0.01;
%     zmin=min(Z(:));zmax=max(Z(:));
%     zmin=dz*floor(zmin/dz);zmax=dz*ceil(zmax/dz);
%     Nmed=nanmedian(histcounts(Z(:),zmin:dz:zmax));
% end
% 
% fprintf('Elevation bin size: %f4.2
% 
% 
% for n=1:size(SG,2)
%     
% xg=SG(n).X;yg=SG(n).Y;zg=SG(n).Z;
% 
% % loop through each survey point to calcualte its z,beta, and xshore
% % deviations
%   for m=5000:5000 %m=1:length(zg)
%       
%       % find the N survey points closest to the same elevation
%       zdiff=abs(zg-zg(m)); % elev diff from point being assessed
%       [zsort,zidx]=sort(zdiff); % sort from closest to furthest
%       zmean=mean(zg(zidx(1:N)));
%       zsigma=std(zg(zidx(1:N)));
%       dsigz(m)=abs(zg(m)-zmean)/zsigma; 
%       
%       fxmean=nanmean(fx(zidx(1:N)));
%       fxsigma=nanstd(fx(zidx(1:N)));
%       dsigfx(m)=abs(fx(m)-fxmean)/fxsigma; 
%           
%       fymean=nanmean(fy(zidx(1:N)));
%       fysigma=nanstd(fy(zidx(1:N)));
%       dsigfy(m)=abs(fy(m)-fymean)/fysigma; 
%       
%       [a,b]=linfit(xg(zidx(1:N)),yg(zidx(1:N)));
%       %d=abs(a.*xg(zidx(1:N)) - yg(zidx(1:N)) + b)./sqrt(a.^2+1);
%       d=(a.*xg(zidx(1:N)) - yg(zidx(1:N)) + b)./sqrt(a.^2+1);
%       dmean=nanmean(d);
%       dsigma=nanstd(d);
%       %dm=abs(a.*xg(m) - yg(m) + b)./sqrt(a.^2+1);
%       dm=(a.*xg(m) - yg(m) + b)./sqrt(a.^2+1);
%       dsigd(m)=abs(dm-dmean)/dsigma;
%   end
% 
% slope=sign(cosd(Mop.Normal(503)-(270-atan2d(fy,fx)))).*sqrt(fx.^2+fy.^2);
% 
% %xg=GM.X2D;yg=GM.Y2D;zg=GM.Z2Dmean;cg=GM.Z2Dclass;
% 
% end
% 
% histogram2(allG,allZ,'facecolor','flat')
% set(gca,'clim',[0 50])
% colormap(jet)
% colorbar