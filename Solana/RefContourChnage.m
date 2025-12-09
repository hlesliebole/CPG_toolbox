% plots shorebox elevation changes all relative to fixed xshore (Y) locations
%  of a specified contour elevation (eg. mhw) on a starting pre- or postnourishment
%  data (eg. 16 apr 2024).

% desired reference contour elevation
%zr=0.774; %msl
%zr=1.344; %mhw
zr=1.566; %mhhw
zr=3.9

% xshore range to keep in rectified data
ymin=-100;
ymax=100;

% load starting survey date X Y Z(ny,nx) matrices 
load SolanaPostnourishmentGrid.mat % 16 apr 2024
[Xr0,Yr0,Zr0,zy]=XshoreRectifyShoreboxBaseline(X,Y,Z,zr,ymin,ymax);

% load ending survey date X Y Z(ny,nx) matrices 
load SolanaPostnourishmentGrid26Apr24.mat
%load SolanaPostnourishmentGrid7Jul24.mat 
[Xr1,Yr1,Zr1,iyzr1]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax);

%figure;surf(Xr,Yr,Zr-zr);polarmap;shading flat;colorbar
figure;surf(Xr1,Yr1,Zr1-Zr0);colormap(flipud(polarmap));shading flat;colorbar
view(2);set(gca,'clim',[-2 2],'ylim',[-20 20],'xlim',[775 2600])
hold on;plot3(Xr1(1,:),iyzr1,iyzr1*0+10,'k.-')

load SolanaPostnourishmentGrid01May24.mat
[Xr2,Yr2,Zr2,iyzr2]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax);

figure;surf(Xr2,Yr2,Zr2-Zr1);colormap(flipud(polarmap));shading flat;colorbar
view(2);set(gca,'clim',[-2 2],'ylim',[-20 20],'xlim',[775 2600])
hold on;plot3(Xr1(1,:),iyzr2,iyzr2*0+10,'k.-')

%%

function [Xr,Yr,Zr,zy]=XshoreRectifyShoreboxBaseline(X,Y,Z,zr,ymin,ymax)

% find xshore location of zr for each shorebox x

for nx=1:size(Z,2)
     iyz=find(Z(:,nx) < zr,1,'first');
    if ~isempty(iyz)
        zy(nx)=iyz;
    else
        zy(nx)=NaN;
    end
end

% make a Yr=ymin:ymax range xshore rectified version of the shorebox X Y Z 
%  where Yr=0 is the reference contour position 

Xr=NaN(numel(ymin:ymax),size(X,2));
Yr=NaN(numel(ymin:ymax),size(Y,2));
Zr=NaN(numel(ymin:ymax),size(Z,2));

for nx=1:size(Z,2)
    if ~isnan(zy(nx))
      Xr(:,nx)=X(zy(nx)+ymin:zy(nx)+ymax,nx);
      Zr(:,nx)=Z(zy(nx)+ymin:zy(nx)+ymax,nx);
      Yr(:,nx)=ymin:ymax;
    end
end

end

%% 

function [Xr,Yr,Zr,iyzr]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax)

% find xshore location of zr for each shorebox x

for nx=1:size(Z,2)
     iyz=find(Z(:,nx) < zr,1,'first');
    if ~isempty(iyz)
        iyzr(nx)=iyz-zy(nx);
    else
        iyzr(nx)=NaN;
    end
end

% make a Yr=ymin:ymax range xshore rectified version of the shorebox X Y Z 
%  where Yr=0 is the reference contour position 

Xr=NaN(numel(ymin:ymax),size(X,2));
Yr=NaN(numel(ymin:ymax),size(Y,2));
Zr=NaN(numel(ymin:ymax),size(Z,2));

for nx=1:size(Z,2)
    if ~isnan(zy(nx))
      Xr(:,nx)=X(zy(nx)+ymin:zy(nx)+ymax,nx);
      Zr(:,nx)=Z(zy(nx)+ymin:zy(nx)+ymax,nx);
      Yr(:,nx)=ymin:ymax;
    end
end

end