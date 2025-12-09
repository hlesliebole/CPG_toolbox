% plots shorebox elevation changes all relative to fixed xshore (Y) locations
%  of a specified contour elevation (eg. mhw) on a starting pre- or postnourishment
%  data (eg. 16 apr 2024).

% desired reference contour elevation
%zr=0.774; %msl
%zr=1.344; %mhw
zr=1.566; %mhhw
zr=3.9;
zr=2.5

% xshore range to keep in rectified data
ymin=-100;
ymax=100;

% shorebox survey files
sf=vertcat(...
'SolanaPostnourishmentGrid26Apr24.mat',...
'SolanaPostnourishmentGrid01May24.mat',...
'SolanaPostnourishmentGrid10May24.mat',...
'SolanaPostnourishmentGrid24May24.mat',...
'SolanaPostnourishmentGrid29May24.mat',...
'SolanaPostnourishmentGrid07Jun24.mat',...
'SolanaPostnourishmentGrid11Jun24.mat',...
'SolanaPostnourishmentGrid21Jun24.mat',...
'SolanaPostnourishmentGrid27Jun24.mat',...
'SolanaPostnourishmentGrid07Jul24.mat');

% load starting survey date X Y Z(ny,nx) matrices 
load SolanaPostnourishmentGrid.mat % 16 apr 2024
[Xr0,Yr0,Zr0,zy]=XshoreRectifyShoreboxBaseline(X,Y,Z,zr,ymin,ymax);

figure('position',[267    63   894   734]);
np=size(sf,1);
nn=0;
for n=1:np
    if n > 1
        Zr1=Zr2;
        iyzr1=iyzr2;
    else
        Zr1=Zr0;
        Xr1=Xr0;
        iyzr1=zy*0;
    end

% load ending survey date X Y Z(ny,nx) matrices 
load(sf(n,:))
%load SolanaPostnourishmentGrid7Jul24.mat 
[Xr2,Yr2,Zr2,iyzr2]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax);

nn=nn+1;
subplot(np,2,nn)
%figure;surf(Xr,Yr,Zr-zr);polarmap;shading flat;colorbar
surf(Xr2,Yr2,Zr2-Zr0);colormap(flipud(polarmap));shading flat;%colorbar
view(2);set(gca,'clim',[-1.5 1.5],'ylim',[-20 20],'xlim',[775 2600])
hold on;
plot3(Xr1(1,:),iyzr2*0,iyzr2*0+10,'w-');
pl=plot3(Xr1(1,:),iyzr2,iyzr2*0+10,'k-');box on;grid on;
if nn < np*2-1;set(gca,'xticklabels',[]);%colorbar('southoutside');
    if n == 1
        text(500,60,'Natural Terrace Edge Elev (2.5m NAVD88 Contour) & Beach Face Change after 16 Apr PostNourishment Survey','fontsize',16);
        title('Cumulative Change','fontsize',16)
        legend(pl,'Natural Terrace Elev','location','south')
    end
    if n == 6
        ylabel('Xshore distance from Natural Terrace Edge Elev Contour on 16 Apr (m)','fontsize',14)
    end
end
if n == np
        xlabel('S to N Alongshore Distance (m)','fontsize',14)
end

nn=nn+1;
subplot(np,2,nn)
%figure;surf(Xr,Yr,Zr-zr);polarmap;shading flat;colorbar
surf(Xr2,Yr2,Zr2-Zr1);colormap(flipud(polarmap));shading flat;%colorbar
view(2);set(gca,'clim',[-1.5 1.5],'ylim',[-20 20],'xlim',[775 2600])
hold on;
plot3(Xr1(1,:),iyzr2*0,iyzr2*0+10,'w-');
plot3(Xr1(1,:),iyzr2,iyzr1*0+10,'-','color',[0 0 0]);box on;grid on;
%plot3(Xr1(1,:),iyzr2-iyzr1,iyzr2*0+10,'k-');box on;
if nn < np*2-1;set(gca,'xticklabels',[]);end
ylabel([sf(n,26:30) ' '],'rotation',0,'fontsize',16)
if n == 1
        title('Incremental Change','fontsize',16)
end
if n == np
        xlabel('S to N Alongshore Distance (m)','fontsize',14)
end

end

cb=colorbar('location','southoutside','fontsize',14);
cb.Position=[0.39    0.05501    0.2729    0.0054];
cb.Label.String='Elevation Change (m)';
%title(cb,'Elevation Change (m)')

% 
% 
% load SolanaPostnourishmentGrid01May24.mat
% [Xr2,Yr2,Zr2,iyzr2]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax);
% 
% figure;surf(Xr2,Yr2,Zr2-Zr1);colormap(flipud(polarmap));shading flat;colorbar
% view(2);set(gca,'clim',[-2 2],'ylim',[-20 20],'xlim',[775 2600])
% hold on;plot3(Xr1(1,:),iyzr2,iyzr2*0+10,'k.-')

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