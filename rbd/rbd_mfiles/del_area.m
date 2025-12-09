
% mfile for deleting data in a defined area      

x=get(gca,'Xlim');
y=get(gca,'Ylim');

jd=find(max(xp) > lon & lon > min(xp) & max(yp) > lat & lat > min(yp));
%kd=find(max(xp) > lon & lon > min(xp));
%jd=find(max(yp) > lat(kd) & lat(kd) > min(yp));
lon(jd)=[];
lat(jd)=[];
dep(jd)=[];
clear kd
clear jd

close(gcf)
make_fig
plot_rbd
