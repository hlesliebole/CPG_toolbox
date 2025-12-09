% buttons to navigate between survey dates

bw=120;db=130;dbs=60;

DN=[SA(ldx).Datenum];
%SA(ldx(CurrentSurveyNumber)).Datenum=SA(ldx(CurrentSurveyNumber)).Datenum;
% [~,CurrentSurveyNumber]=min(abs(DN-(SA(ldx(CurrentSurveyNumber)).Datenum+365)));
% datestr(SA(ldx(CurrentSurveyNumber)).Datenum)
% datestr(SA(ldx(CurrentSurveyNumber)).Datenum)

%%
% SurveyDatetime=datetime([SA.Datenum],'ConvertFrom','datenum');
% CurrentYear=year(SurveyDatetime(CurrentSurveyNumber));
% CurrentMonth=month(SurveyDatetime(CurrentSurveyNumber));
% PrevYearSurveyNumber=find(year(SurveyDatetime) == CurrentYear-1,1,'first');
% if isempty(PrevYearSurveyNumber);PrevYearSurveyNumber=1;end
% NextYearSurveyNumber=find(year(SurveyDatetime) == CurrentYear+1,1,'first');
% if isempty(NextYearSurveyNumber);NextYearSurveyNumber=size(SA,2);end

% MainFig=figure('position',[250  103  1183 694]);
% ax1=axes('position',[.1 .1 .8 .8]);

nav(1)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs 655 bw 35],...
'Callback','CurrentSurveyNumber=1;PlotSAQC',...
'String','First Surv'); set(nav(1),'ForegroundColor','blue');

nav(2)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db 655 bw 35],...
'Callback','[~,CurrentSurveyNumber]=min(abs(DN-(SA(ldx(CurrentSurveyNumber)).Datenum-365)));PlotSAQC',...
'String','<- Prev Year'); set(nav(2),'ForegroundColor','blue');

nav(3)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*2 655 bw 35],...
'Callback','[~,CurrentSurveyNumber]=min(abs(DN-(SA(ldx(CurrentSurveyNumber)).Datenum-30.5)));PlotSAQC',...
'String','<- Prev Month'); set(nav(3),'ForegroundColor','blue');

nav(4)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*3 655 bw 35],...
'Callback','CurrentSurveyNumber=max([1 CurrentSurveyNumber-1]);PlotSAQC',...
'String','<- Prev Surv'); set(nav(4),'ForegroundColor','blue');

nav(5)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*4 655 bw 35],...
'Callback','CurrentSurveyNumber=min([numel(DN) CurrentSurveyNumber+1]);PlotSAQC',...
'String','Next Surv ->'); set(nav(5),'ForegroundColor','red');

nav(6)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*5 655 bw 35],...
'Callback','[~,CurrentSurveyNumber]=min(abs(DN-(SA(ldx(CurrentSurveyNumber)).Datenum+30.5)));PlotSAQC',...
'String','Next Month ->'); set(nav(6),'ForegroundColor','red');

nav(7)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*6 655 bw 35],...
'Callback','[~,CurrentSurveyNumber]=min(abs(DN-(SA(ldx(CurrentSurveyNumber)).Datenum+365)));PlotSAQC',...
'String','Next Year ->'); set(nav(7),'ForegroundColor','red');

nav(8)=uicontrol(gcf,'Style','Pushbutton','Position',[dbs+db*7 655 bw 35],...
'Callback','CurrentSurveyNumber=numel(DN);PlotSAQC',...
'String','Last Surv'); set(nav(8),'ForegroundColor','red');

