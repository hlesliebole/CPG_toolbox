
% mfile for selecting a bathymetry point to add


zoom off

htxt=uicontrol(gcf,'Style','text','Position',[10 top-30 100 20],...
'String','Depth'); 

hed=uicontrol(gcf,'Style','edit','Position',[10 top-60 100 30],...
'String','0'); 

% make dummy button for quitting

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','igo=1',...
'String','Done'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 10 500 20],...
'String','Using cursor, edit depth value and then select location '); 

for idum=1:5000
xy=ginput(1);

v=get(gca,'Ylim');
if xy(2) < v(2)       
if xy(2) > v(1)       

depth=str2num(get(hed,'String'))
if -depth > dmax, depth=-(dmax-.1), end
if(depth >= 0);
nc=1;
else
nc=1+fix(31.*log(-depth)/log(dmax));
end
%nc=1+fix(32.*log(-depth)/log(dmax));
if nc < 1, nc=1 , end;

' new point list is now : '
nn=nn+1;
new(nn,1)=xy(1);
new(nn,2)=xy(2);
new(nn,3)=depth;
new

% add new point to lat lon dep
ne=max(size(lon));
ne=ne+1;

lon(ne)=xy(1);
lat(ne)=xy(2);
dep(ne)=depth;


% plot new point

if nn > 0;
 plot(new(nn,1),new(nn,2),'.','Color',c32(nc,1:3),'MarkerSize',ms)
end;

end;
end;

if xy(2) < v(1)
zoom on
delete(hpop)
delete(htxt)
delete(hed)
%close(gcf)
%make_fig
%plot_rbd
return
end

end

