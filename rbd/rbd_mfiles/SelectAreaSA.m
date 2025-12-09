
% if in 3d view, switch to top mode to find points with ginput
if dview == '3'
  dview = 'T';
  PlotSAQC
end

if exist('nav','var');delete(nav);end
% mfile for setting area to delete 
zoom on 

hpop=uicontrol(gcf,'Style','Pushbutton','Position',[10 625 100 50],...
'Callback','DeleteAreaSA','String','Delete'); 
set(hpop,'ForegroundColor','red');

htxt=uicontrol(gcf,'Style','text','Position',...
[170 670 800 25],'String', ...
'Zoom in (click or click-drag mouse) to define delete area. Then select Delete button.'); 
set(htxt,'ForegroundColor','red','fontsize',18);

hpop2=uicontrol(gcf,'Style','Pushbutton','Position',[10 10 100 50],...
'Callback','delete(hpop);delete(htxt);delete(hpop2);PlotSAQC;edit_menuSA','String','Done'); 

% 
% 
% % mfile for deleting data in a defined area      
% 
% xla=get(gca,'Xlim');
% yla=get(gca,'Ylim');
% 
% I=find(x < max(xla)  & x > min(xla) & y < max(yla) & y > min(yla));
% %kd=find(max(xp) > lon & lon > min(xp));
% %jd=find(max(yp) > lat(kd) & lat(kd) > min(yp));
% nb=numel(QC(CurrentSurveyNumber).X);
% nadd=nb+1:nb+numel(I);
% QC(CurrentSurveyNumber).X(nadd)=x(I);
% QC(CurrentSurveyNumber).Y(nadd)=y(I);
% QC(CurrentSurveyNumber).Z(nadd)=z(I);
% 
% % remove duplicate deleted points
% % check for duplicate x, y, z
% [uxy,ia,ic]=unique([[QC(CurrentSurveyNumber).X]',[QC(CurrentSurveyNumber).Y]'],'rows');
% % keep unique points
% QC(CurrentSurveyNumber).X=QC(CurrentSurveyNumber).X(ia);
% QC(CurrentSurveyNumber).Y=QC(CurrentSurveyNumber).Y(ia);
% QC(CurrentSurveyNumber).Z=QC(CurrentSurveyNumber).Z(ia);
% 
% % mark point as deleted
% 
% %if nb > 0
%     hold on;
%     if exist('dp','var');delete(dp);end
%  if dview == '3'
%     dp=plot3([QC(CurrentSurveyNumber).X],[QC(CurrentSurveyNumber).Y],...
%      [QC(CurrentSurveyNumber).Z],'mx','MarkerSize',10,'linewidth',2);
% elseif dview == 'T'
%     dp=plot([QC(CurrentSurveyNumber).X],[QC(CurrentSurveyNumber).Y],...
%      'mx','MarkerSize',10,'linewidth',2);
% elseif dview == 'N'
%     dp=plot([QC(CurrentSurveyNumber).X],[QC(CurrentSurveyNumber).Z],...
%      'mx','MarkerSize',10,'linewidth',2);
% elseif dview == 'E'
%     dp=plot([QC(CurrentSurveyNumber).Y],[QC(CurrentSurveyNumber).Z],...
%      'mx','MarkerSize',10,'linewidth',2);
% end
% %end
% 
% % loop around
% 
% delete(hpop)
% delete(htxt)
% DeleteAreaSA
