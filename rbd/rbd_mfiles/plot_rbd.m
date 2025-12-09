%
% Makes a color dot plot of nos bathymetry data read in by
% load_rbd.m
%
%

% plot the coastline if you have one

if exist('coast.xy')
 load coast.xy
 plot(coast(1,:),coast(2,:),'w-')
end

% set max depth
min_depth=min(-dep)
max_depth=max(-dep)
dmax=max(-dep)+1;

% plot the soundings color coded by depth
%  color scale selected in defaults.m

c32=c(1:2:64,:);

' Sorting depths into 32 colors '

x=1+31*log(-dep)/log(dmax);
x(1)
% fix log of 0 depths
x(isinf(x) == 1)=1;

%for n=32:-1:1
for n=1:1:32
%xx=n-floor(x);
%xx=sign(xx);
%x=abs(sign(x));
%x=abs(x-1);
y=find(x >= n & x <= (n+1));
col=c32(n,1:3);
plot(lon(y),lat(y),'.','Color',col,'MarkerSize',ms)
end

if (exist('spot.dat','file') == 2 );
fid=fopen('spot.dat','r');
spll=fscanf(fid,'%g %g\n');
fclose(fid)
plot(spll(1),spll(2),'m*','markersize',10);
plot(spll(1),spll(2),'ko','markersize',10);
end

clear x;
clear y;

set(gcf,'Visible','on')
