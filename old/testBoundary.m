%close all
load M00503SA.mat
n=6;
x=SA(n).X;y=SA(n).Y;z=SA(n).Z;

%----------------------------------------------------------
% Step 1: Put x,y,z spatially averaged data in a 2D elevation 
%         grid arrays, X,Y,Z. "no data" grid points = NaNs.
%----------------------------------------------------------

nx=max(x)-min(x)+3;ny=max(y)-min(y)+3; % grid dimensions
Z=nan(ny,nx); % initial grid of NaNs
gidx=sub2ind([ny nx],y-min(y)+2,x-min(x)+2); % data point 1d grid indices
Z(gidx)=z; % add elevation data to temp grid
Z=round(Z,3);
[X,Y]=meshgrid(min(x)-1:max(x)+1,min(y)-1:max(y)+1); % 2D grid indices
clf;subplot(1,2,1);surf(X,Y,Z);colormap(jet);set(gca,'clim',[0 4]);%demcmap(Z);
xl=get(gca,'xlim');yl=get(gca,'ylim');zl=get(gca,'zlim');

%[np,edges]=histcounts(Z(:),[.01*floor(100*min(Z(:))):0.01:max(Z(:))]);

%----------------------------------------------------------
% Step 2: Replace any NaNs within Z grid array's survey spatial 
%         bounds with interpolated elevations.
%----------------------------------------------------------

iData=find(~isnan(Z(:))); % 1d indices of grid points with data 
iNoData=find(isnan(Z(:))); % 1d indices of grid points without data
% delaunay interpolate to fill in any no data gaps in survey area
Z(iNoData)=griddata(X(iData),Y(iData),Z(iData),X(iNoData),Y(iNoData));

%-------------------------------------------------------------
% Step 3: Define extreme quantile exceedance differences DqZ
%         at all grid points.
%-------------------------------------------------------------

DpZ=nan(ny,nx);
% get quantile function and min-max slope break probabilities
[qZ,p,pminz,pmaxz]=GetQuantCutoffs(Z,1); 
[qZu,ia,ic]=unique(qZ);pu=p(ia); % remove duplicate function values
% interpolate quantile function p,qZ
pZ=nan(ny,nx);pZ(:)=interp1(qZu,pu,Z(:));
% probability differences from min-max slope break points
DpZ=pZ;DpZ(:)=max([pminz-pZ(:)';pZ(:)'-pmaxz]);
Z2=nan(ny,nx); 
Z2(gidx)=Z(gidx);Z2(DpZ > 0)=NaN;
%figure;surf(Z2);colormap(jet);set(gca,'clim',[0 4]);%demcmap(Z);


Nmin=500; % min number of points allowed to define a elevation boundary area
%dp1=0.01;
dp=0.01;
% p1=0-dp1;
p2=0;
%p1=0-dp1;
p2=dp/8;


Zsave=Z;
%zmax=0;
 while p2 < 1
    %p1=p1+dp1;p2=p1;np=0;
    %p1=p2;np=0;
    p1=p2-dp/8;np=0;
    while np < Nmin %&& p2 < 1
        p2=p2+dp;
        % move p1 back if p2 has reached 1 and there are too few points
        if(p2 >= 1);p2=1;p1=p1-dp;end 
        zmin=quantile(Z(:),p1);zmax=quantile(Z(:),p2);
        idx=find(Z(:) >= zmin & Z(:) <= zmax);
        
        np=length(idx);
    end
    
     k=boundary(X(idx),Y(idx));nk=length(k);
%        
%    % get quantile function and min-max slope break probabilities
%    [qZ,p,pminz,pmaxz]=GetQuantCutoffs(Z(idx(k)),1);
%    [qZu,ia,ic]=unique(qZ);pu=p(ia); % remove duplicate function values
%     % interpolate quantile function p,qZ
%         pZ(idx(k))=interp1(qZu,pu,Z(idx(k)));
%         %pZ=nan(ny,nx);pZ(:)=interp1(qZ,p,Z(zip));
%    % probability differences from min-max slope break points
%    %DpZ=pZ;
%    DpZ(idx(k))=max([pminz-pZ(idx(k))';pZ(idx(k))'-pmaxz]);
%    nE=length(find(DpZ(idx(k)) > 0));
%    Z(DpZ > 0)=NaN;
%    nE=length(find(Z(idx(k)) > zmax));
%      
%      
%      
     in=inpolygon(X,Y,X(idx(k)),Y(idx(k)));
     zip=find(in(:) == 1);nb=length(zip);
     
     figure;
     histogram(Z(zip),(.01*floor(100*min(Z(zip)))-0.01):0.01:(max(Z(zip))+0.01));
     [np,edges]=histcounts(Z(zip),(.01*floor(100*min(Z(zip)))-0.01):0.01:(max(Z(zip))+0.01));
     % find zero bins closest to the peak bin in each direction
     [bmax,npk]=max(np);nzero=find(np == 0);ndist=nzero-npk;
     dmin=max(ndist(ndist < 0));zbmin=edges(nzero(ndist == dmin)+1);
     dmax=min(ndist(ndist > 0));zbmax=edges(nzero(ndist == dmax));
     
     
     if(nanmean(Z(zip)) < zmin)
            zmin=nanmean(Z(zip));
            fprintf('Adjusting zmin down: %6.3f\n',zmin)
            idx=find(Z(:) >= zmin & Z(:) <= zmax);
            in=inpolygon(X,Y,X(idx(k)),Y(idx(k)));
            zip=find(in(:) == 1);nb=length(zip);
        end
        if(nanmean(Z(zip)) > zmax)   
            zmax=nanmean(Z(zip));
            fprintf('Adjusting zmax up: %6.3f\n',zmax)
            idx=find(Z(:) >= zmin & Z(:) <= zmax);
            in=inpolygon(X,Y,X(idx(k)),Y(idx(k)));
            zip=find(in(:) == 1);nb=length(zip);
        end
      
     nE=length(find(Z(zip)-zmax > zmax-zmin));
     Emax=max(Z(zip)-zmax);
     Emin=min(Z(zip)-zmin);
     %zip=unique([zip' idx(k)']);
     
     fprintf('%5.3f %5.3f %6.3f %6.3f %6i %6i %6i %6i %8.3f %8.3f\n',...
         p1,p2,zmin,zmax,np,nk,nb,nE,Emin,Emax)
     
     %fprintf('sigma/dz ratio %6.3f\n',nanstd(Z(zip))/(zmax-zmin))
     
     if nE > 0 && p2 < 1 
         
   % get quantile function and min-max slope break probabilities
   [qZ,p,pminz,pmaxz]=GetQuantCutoffs(Z(zip),1);
   %pminz=max([pminz 0.0025]);pmaxz=min([pmax 0.9975]);
   %pminz=min([pminz 0.01]);pmaxz=max([pmax 0.99]);
   [qZu,ia,ic]=unique(qZ);pu=p(ia); % remove duplicate function values
    % interpolate quantile function p,qZ
        pZ=interp1(qZu,pu,Z(zip));
        pmin4sigma=interp1(qZu,pu,nanmean(Z(zip))-4*nanstd(Z(zip)));
        pmax4sigma=interp1(qZu,pu,nanmean(Z(zip))+4*nanstd(Z(zip)));
        
        if pminz > pmin4sigma;fprintf('Using Quant min\n');else
            fprintf('Using 4Sigma min  %6.3f %6.3f %6.3f %6.3f %6.3f\n',...
                min(Z(zip)),max(Z(zip)),nanmean(Z(zip)),nanstd(Z(zip)),nanmean(Z(zip))-4*nanstd(Z(zip)));end
        if pmaxz < pmin4sigma;fprintf('Using Quant max\n');else
            fprintf('Using 4Sigma max %6.3f %6.3f %6.3f %6.3f %6.3f\n',...
                min(Z(zip)),max(Z(zip)),nanmean(Z(zip)),nanstd(Z(zip)),nanmean(Z(zip))+4*nanstd(Z(zip)));end
        pminz=max([pminz pmin4sigma]);pmaxz=min([pmax pmax4sigma]);
        
        fprintf('3sigma %6.3f %6.3f\n',nanmean(Z(zip))+3*nanstd(Z(zip)),max(Z(zip)))
        
        if max(Z(zip)) <= nanmean(Z(zip))+3*nanstd(Z(zip))
            pmaxz=1;
        end
        
        if min(Z(zip)) >= nanmean(Z(zip))-3*nanstd(Z(zip))
            pminz=0;
        end
        
        
        
        
        %quantile(Z(zip),pmaxz)
        %pZ=nan(ny,nx);pZ(:)=interp1(qZ,p,Z(zip));
   % probability differences from min-max slope break points
   %DpZ=pZ;
   DpZtmp=max([pminz-pZ';pZ'-pmaxz]);
   
   fprintf('%i points flagged\n',length(find(DpZtmp > 0)));
   DpZ(zip)=nanmax([DpZ(zip)';DpZtmp]);
   Z(DpZ > 0)=NaN;
   
iData=find(~isnan(Z(:))); % 1d indices of grid points with data 
iNoData=find(isnan(Z(:))); % 1d indices of grid points without data
% delaunay interpolate to fill in any no data gaps in survey area
Z(iNoData)=griddata(X(iData),Y(iData),Z(iData),X(iNoData),Y(iNoData));
     end
 
 end

%figure;surf(X,Y,Z);colormap(jet);set(gca,'clim',[0 4]);
 
Z=Zsave;
Z2=nan(ny,nx); 
Z2(gidx)=Z(gidx);Z2(DpZ > 0)=NaN;subplot(1,2,2);surf(X,Y,Z2);colormap(jet);set(gca,'clim',[0 4]);%demcmap(Z);
set(gca,'xlim',xl);set(gca,'ylim',yl);set(gca,'zlim',zl);
% figure;surf(X,Y,Z);view(3);demcmap(Z);hold on;
% hold on;plot3(X(idx(k)),Y(idx(k)),Z(idx(k))+.05,'g-','linewidth',2);
% d=sqrt(diff(X(idx(k))).^2+diff(Y(idx(k))).^2);
% 
% j=find(in(:) == 1);length(j)
% %plot3(X(j),Y(j),Z(j)+0.1,'c.','markersize',5);
% l=find(Z(j) < zmin | Z(j) > zmax );
% plot3(X(j(l)),Y(j(l)),Z(j(l)),'r.','markersize',10);
% % figure;
% % plot3(X(idx(j)),Y(idx(j)),'k.','markersize',3);