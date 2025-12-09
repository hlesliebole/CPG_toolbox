% delete points in the area shown on the plot axes

delete(hpop);
delete(htxt);
delete(hpop2);

dlx=[QC(ldx(CurrentSurveyNumber)).X];
dly=[QC(ldx(CurrentSurveyNumber)).Y];
dlz=[QC(ldx(CurrentSurveyNumber)).Z];

xla=get(gca,'Xlim');
yla=get(gca,'Ylim');

if dview == 'T'
  I=find(dlx < max(xla)  & dlx > min(xla) & dly < max(yla) & dly > min(yla));
end
% looking north, match with x and z
if dview == 'N'
 I=find(dlx < max(xla) & dlx > min(xla) & dlz < max(yla) & dlz > min(yla));
end
% looking east, match with y and z
if dview == 'E'
 I=find(dly < max(xla) & dly > min(xla) & dlz < max(yla) & dlz > min(yla));
end

% remove point from delete list
QC(ldx(CurrentSurveyNumber)).X(I)=[];
QC(ldx(CurrentSurveyNumber)).Y(I)=[];
QC(ldx(CurrentSurveyNumber)).Z(I)=[];

% replot the data
PlotSAQC

% mark all the delated points
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

UNSelectAreaSA
