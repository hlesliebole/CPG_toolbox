load SolanaPrenourishmentGrid.mat
load SolanaPostnourishmentGrid.mat
load SolanaShoreboxMap.mat
idx=find(Z(:) > -12.5 & ~isnan(Z0(:)));
dZ=Z*NaN;
dZ(idx)=Z(idx)-Z0(idx);
sum(Z(idx)-Z0(idx))
figure
%imagesc(xgimg(1050:3000),ygimg,Z(:,1050:3000));
imagesc(xgimg,ygimg,dZ,'alphadata',(~isnan(dZ)));
colorbar
set(gca,'clim',[-5 5],'ydir','normal')
colormap(flipud(polarmap));

figure;plot(sum(dZ,'omitnan'))

iz0=round(Z0(idx)*10);
dz0=dZ(idx);
[siz0,idz0]=sort(iz0);
sdz0=dZ(idz0);

%%
% vertical resolution for summed changes
Res=0.1;
% round survey x,y to desired resolution
%zr=Res*round(z1da/Res); % round to Res meter spatial resolution
zr=Res*round(Z0(:)/Res);

dz=dZ(:);

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');
