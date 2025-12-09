
%close(gcf)

%sort_rbd
%dup_rbd
%write_rbd

%!tgrid_rbd

[xq,yq]=meshgrid(xp(1):0.001:xp(2),yp(1):0.001:yp(3));
vq=griddata(lon,lat,dep,xq,yq);
pcolor(xq,yq,vq);hold on;shading interp;
%p=plot(lon,lat,'k.');uistack(p,'top')
set(gca,'clim',[-100 0]);colorbar;
