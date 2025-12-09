
% mfile for setting area to delete 

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','show_grid','String','Delete'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 10 500 20],...
'String','Define area by drag and drop with left mouse button'); 

% get rid of aspect ratio 

%set(h,'DataAspectRatio',[NaN, NaN])

zoom on
