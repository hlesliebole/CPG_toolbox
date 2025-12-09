
% mfile for selecting a bathymetry point to delete

zoom off

% make dummy button for quitting

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','delete(hpop)','String','Done'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 10 500 20],...
'String','Select points to delete with the cursor'); 

xy=ginput(1);

v=get(gca,'Ylim')
if xy(2) < v(1)       
zoom on
return
end

% find closet lat lon point
z=abs(lon-xy(1))+abs(lat-xy(2));
[Y,I]=min(z);


' Bad list is now : '
nb=nb+1;
bad(nb,1)=lon(I);
bad(nb,2)=lat(I);
bad(nb,3)=dep(I);
bad

% remove data

lon(I)=[];
lat(I)=[];
dep(I)=[];

% mark point as deleted

if nb > 0;
 plot(bad(:,1),bad(:,2),'w*','MarkerSize',ms)
end;

% loop around

delete(hpop)
del_rbd
