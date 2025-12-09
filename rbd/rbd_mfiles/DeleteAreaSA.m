
%% mfile for deleting data in a defined area  

delete(hpop);
delete(htxt);
delete(hpop2);

xla=get(gca,'Xlim');
yla=get(gca,'Ylim');

if dview == 'T'
  I=find(x < max(xla)  & x > min(xla) & y < max(yla) & y > min(yla));
end
% looking north, match with x and z
if dview == 'N'
 I=find(x < max(xla) & x > min(xla) & z < max(yla) & z > min(yla));
end
% looking east, match with y and z
if dview == 'E'
 I=find(y < max(xla) & y > min(xla) & z < max(yla) & z > min(yla));
end

nb=numel(QC(ldx(CurrentSurveyNumber)).X);
nadd=nb+1:nb+numel(I);
QC(ldx(CurrentSurveyNumber)).X(nadd)=x(I);
QC(ldx(CurrentSurveyNumber)).Y(nadd)=y(I);
QC(ldx(CurrentSurveyNumber)).Z(nadd)=z(I);

% remove any duplicate deleted points
% check for duplicate x, y, z
[uxy,ia,ic]=unique([[QC(ldx(CurrentSurveyNumber)).X]',[QC(ldx(CurrentSurveyNumber)).Y]'],'rows');
% keep unique points
QC(ldx(CurrentSurveyNumber)).X=QC(ldx(CurrentSurveyNumber)).X(ia);
QC(ldx(CurrentSurveyNumber)).Y=QC(ldx(CurrentSurveyNumber)).Y(ia);
QC(ldx(CurrentSurveyNumber)).Z=QC(ldx(CurrentSurveyNumber)).Z(ia);
    
% mark points as deleted

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
%end

SelectAreaSA
% loop around

% delete(hpop)
% delete(htxt)
% DeleteAreaSA
