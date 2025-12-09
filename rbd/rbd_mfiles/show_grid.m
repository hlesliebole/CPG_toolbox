
% mfile for running the gridding program 

x=get(gca,'Xlim');
y=get(gca,'Ylim');

x=-x;
x
y
save tgrid_rbd.inp x y -ascii
' Grid limits saved in tgrid_rbd.inp'

delete(hpop)
delete(htxt)
set(h,'DataAspectRatio',[rat 1])
zoom out
x=-x;
xp=[x(1) x(2) x(2) x(1) x(1)];
yp=[y(1) y(1) y(2) y(2) y(1)];
plot(xp,yp,'k-');

