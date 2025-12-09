
% if in 3d view, switch to top mode to find points with ginput
if dview == '3'
  dview = 'T';
  PlotSAQC
end

% remove survey nav buttons
if exist('nav','var');delete(nav);end

% mfile for selecting a bathymetry point to delete

zoom off

% make dummy button for quitting

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','delete(hpop);delete(htxt);PlotSAQC;edit_menuSA','String','Done'); 

htxt=uicontrol(gcf,'Style','text','Position',...
[150 670 800 25],'String', ...
'Click on individual points to delete them. Click on Done button when finished.'); 
set(htxt,'ForegroundColor','red','fontsize',18);

xy=ginput(1);

v=get(gca,'Ylim');
if xy(2) < v(1)       
delete(hpop)
delete(htxt)
PlotSAQC;
edit_menuSA
return
end

% find closet lat lon point
%vw=get(ax1,'view');
% looking from top, match with x and y
if dview == 'T'
 d=sqrt((x-xy(1)).^2+(y-xy(2)).^2);
 [Y,I]=min(d);
end
% looking north, match with x and z
if dview == 'N'
 d=sqrt((x-xy(1)).^2+(z-xy(2)).^2);
 [Y,I]=min(d);
end
% looking east, match with y and z
if dview == 'E'
 d=sqrt((y-xy(1)).^2+(z-xy(2)).^2);
 [Y,I]=min(d);
end


%' Bad list is now : '
nb=numel(QC(ldx(CurrentSurveyNumber)).X);
nb=nb+1;
QC(ldx(CurrentSurveyNumber)).X(nb)=x(I);
QC(ldx(CurrentSurveyNumber)).Y(nb)=y(I);
QC(ldx(CurrentSurveyNumber)).Z(nb)=z(I);

% remove duplicate deleted points
% check for duplicate x, y, z
[uxy,ia,ic]=unique([[QC(ldx(CurrentSurveyNumber)).X]',[QC(ldx(CurrentSurveyNumber)).Y]'],'rows');
% keep unique points
QC(ldx(CurrentSurveyNumber)).X=QC(ldx(CurrentSurveyNumber)).X(ia);
QC(ldx(CurrentSurveyNumber)).Y=QC(ldx(CurrentSurveyNumber)).Y(ia);
QC(ldx(CurrentSurveyNumber)).Z=QC(ldx(CurrentSurveyNumber)).Z(ia);
    
% mark point as deleted

%if nb > 0
    hold on;
    if exist('dp','var');delete(dp);end
 if dview == '3'
    dp=plot3([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Y],...
     [QC(ldx(CurrentSurveyNumber)).Z],'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'T'
    dp=plot([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Y],...
     'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'N'
    dp=plot([QC(ldx(CurrentSurveyNumber)).X],[QC(ldx(CurrentSurveyNumber)).Z],...
     'mx','MarkerSize',10,'linewidth',2);
elseif dview == 'E'
    dp=plot([QC(ldx(CurrentSurveyNumber)).Y],[QC(ldx(CurrentSurveyNumber)).Z],...
     'mx','MarkerSize',10,'linewidth',2);
end
%end

% loop around

delete(hpop)
delete(htxt)
DeletePointsSA
