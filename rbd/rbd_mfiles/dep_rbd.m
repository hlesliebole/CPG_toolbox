
% mfile for showing the depth of a point


zoom off

htxt=uicontrol(gcf,'Style','text','Position',[10 top-30 100 20],...
'String','Depth'); 

hed=uicontrol(gcf,'Style','text','Position',[10 top-60 100 30],...
'String','?'); 

% make dummy button for quitting

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','igo=1',...
'String','Done'); 


for idum=1:5000
xy=ginput(1);

v=get(gca,'Ylim');
if xy(2) < v(2)       
if xy(2) > v(1)       

% find closet lat lon point
z=abs(lon-xy(1))+abs(lat-xy(2));
[Y,I]=min(z);
depth=dep(I) 

delete(hed)
hed=uicontrol(gcf,'Style','text','Position',[10 top-60 100 30],...
'String',num2str(depth)); 

end;
end;

if xy(2) < v(1)
zoom on
delete(hpop)
delete(htxt)
delete(hed)
return
end

end

