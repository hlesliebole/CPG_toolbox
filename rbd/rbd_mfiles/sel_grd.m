
% mfile for running the gridding program 

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','show_grid','String','Save'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 10 500 20],...
'String','Define area by drag and drop with left mouse button'); 

% get rid off aspect ratio 

set(h,'AspectRatio',[NaN NaN])

zoom on
