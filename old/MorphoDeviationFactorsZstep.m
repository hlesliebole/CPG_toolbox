Mopnum=503;

load M00503SA.mat
n=6;
xg=SA(n).X;yg=SA(n).Y;zg=SA(n).Z;
% look at surface gradients of a gridded survey

nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
Z=nan(ny,nx); % temp grid of Nans

idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
Z(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);
% %[fx,fy]=gradient(Z);fx=fx(idx);fy=fy(idx);

igood=find(~isnan(Z(:)));ibad=find(isnan(Z(:)));
Z(ibad)=griddata(X(igood),Y(igood),Z(igood),X(ibad),Y(ibad));

for loop=1:1

% dx1=diff(Z,1,2);dx1(:,end+1)=dx1(:,end); % x grid elev differences
% dy1=diff(Z,1,1);dy1(end+1,:)=dy1(end,:); % y grid elevation differences
% sign1=sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy1,dx1))));
% dx2=fliplr(diff(fliplr(Z),1,2));dx2(:,end+1)=dx2(:,end);
% dy2=flipud(diff(flipud(Z),1,2));dy2(:,end+1)=dy2(:,end);
% sign2=-sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy2,dx2))));

dx1=diff(Z,1,2);dx1(:,end+1)=dx1(:,end); % x grid elev differences
dy1=diff(Z,1,1);dy1(end+1,:)=dy1(end,:); % y grid elevation differences
sign1=sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy1,dx1))));
dx2=diff(fliplr(Z),1,2);dx2(:,end+1)=dx2(:,end);dx2=fliplr(dx2);
dy2=diff(flipud(Z),1,1);dy2(end+1,:)=dy2(end,:);dy2=flipud(dy2);
sign2=-sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy2,dx2))));

signc=sign(sign1.*sqrt(dx1.^2+dy1.^2)+sign2.*sqrt(dx2.^2+dy2.^2));
stpx=(abs(dx1)+abs(dx2))/2;
stpy=(abs(dy1)+abs(dy2))/2;
stp=signc.*sqrt(stpx.^2+stpy.^2);


% elevation step size
dz=0.1;
% minimum number of data points to make qc distribution
Nmin=100;

% figure out step range that captures all the data
minz=dz*floor(min(Z(:)/dz));
maxz=dz*floor(max(Z(:)/dz));

for z1=minz:dz:maxz
    
    z2=z1+dz;
    ielv=find(Z(:) >= z1 & Z(:) < z2);  
    
   if( length(ielv) < Nmin) % fix last sample bin
      zdiff=abs(Z(:)-mean([z1 z2])); % elev diff from point being assessed
      [zsort,zidx]=sort(zdiff); % sort from closest to furthest
      idx=zidx(1:Nmin); % highest N elevation points
      Nmed=Nmin;
      %fprintf('%5.3f %5.3f %10i Expanded Elevation Search\n',...
      %    min(Z(idx)),max(Z(idx)),Nmed)
   else
      idx=ielv;
      Nmed=length(idx);
      fprintf('%5.3f %5.3f %10i\n',z1,z2,Nmed)
   end
   StpDev(idx)=stp(idx); 
   
    % Mean Steepness Deviation
    stpmean=nanmean(stp(idx)); % mean steepness of elev range
    stpsigma=nanstd(stp(idx)); % steepness standard dev sigma of elev range
    StpDev(idx)=abs(stp(idx)-stpmean)/stpsigma; % steepness deviations in fraction of sigma
%       
%     [a,b]=linfit(X(idx),Y(idx)); % best fit line to elv points
%     p1=plot(xl,xl*a+b,'k-');hold on;
%     d=(a.*X(idx) - Y(idx) + b)./sqrt(a.^2+1); % distance from line
%     dmean=nanmean(d); % mean distance from best fit line
%     dsigma=nanstd(d); % standard dev sigma distance
%     XshDev(idx)=(d-dmean)/dsigma;
      
    end
% 
% Na=[Na N];
% S2=[S2 length(find(StpDev > 2))];
% S3=[S3 length(find(StpDev > 3))];
% S4=[S4 length(find(StpDev > 4))];
% 
% 
% fprintf('%5i %10i %10i %10i\n',N,length(find(StpDev > 2)),...
%    length(find(StpDev > 3)),length(find(StpDev > 4))) 
% end

% figure;pcolor(abs(StpDev));shading flat;colormap(jet);colorbar
% set(gca,'clim',[0 3.5]);colormap(jet(7));

Z2=Z;
igood=find(StpDev(:) < 1.5);
istp=find(StpDev(:) >= 1.5);
ibad=istp;%ibad=find(StpDev(:) >= 3 | isnan(StpDev(:)));
fprintf('%i points > 3-sigma steepness\n',length(istp))
Z2(ibad)=nan;
% figure;ScatterPlotMopUTM(X(:),Y(:),Z2(:),'2d')
% S2=StpDev;S2(ibad)=nan;
% figure;ScatterPlotMopUTM(X(:),Y(:),S2(:),'2d')
% set(gca,'color',[.8 .8 .8]);set(gca,'clim',[0 3.5]);colormap(jet(7));

Z2(ibad)=griddata(X(igood),Y(igood),Z2(igood),X(ibad),Y(ibad));
Z=Z2;
end


figure;ScatterPlotMopUTM(X(:),Y(:),Z2(:),'3d');set(gca,'clim',[0 3.5]);

figure;pcolor(abs(StpDev));shading flat;colormap(jet);colorbar
set(gca,'clim',[0 4.5]);colormap(jet(9));


% 
% figure;pcolor(abs(XshDev));shading flat;colormap(jet);colorbar
% set(gca,'clim',[0 3.5]);colormap(jet(7));
% 
% % poor mans fix to truncated normal distribution
% figure; histogram([stp(:)' (2*nanmean(stp(:))-stp(ineg))'])


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