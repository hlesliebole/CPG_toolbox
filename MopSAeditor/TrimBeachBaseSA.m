
%% mfile for deleting edge data with positive west or south slopes

zmax=quantile(z,.05); % only delete points at low end of z elevations
BaseDepth=5;

nb=numel(QC(ldx(CurrentSurveyNumber)).X);

% flag low beach bad slope points going west to east

for yd=min(y):max(y)
    ndx=find(y == yd);
    for nys=ndx
    for xd=min(x(nys)):min(x(nys))+BaseDepth-1
        mdx1=find( y == yd & x == xd );
        mdx2=find( y == yd & x == xd + 1 );
        if z(mdx1) < zmax & ( isempty(mdx2) | z(mdx1) > z(mdx2) )  
            I=mdx1;
            nb=nb+1;
            QC(ldx(CurrentSurveyNumber)).X(nb)=x(I);
            QC(ldx(CurrentSurveyNumber)).Y(nb)=y(I);
            QC(ldx(CurrentSurveyNumber)).Z(nb)=z(I);
        end
    end
    end
end

% 
% % flag low beach bad slope points going south to north
% 
% for xd=min(x):max(x)
%     ndx=find(x == xd);
%     for nys=ndx
%     for yd=min(y(nys)):min(y(nys))+BaseDepth-1
%         mdy1=find( y == yd & x == xd );
%         mdy2=find( y == yd + 1 & x == xd );
%         if z(mdy1) < 1 & ( isempty(mdy2) | z(mdy1) > z(mdy2) )  
%             I=mdy1;
%             nb=nb+1;
%             QC(ldx(CurrentSurveyNumber)).X(nb)=x(I);
%             QC(ldx(CurrentSurveyNumber)).Y(nb)=y(I);
%             QC(ldx(CurrentSurveyNumber)).Z(nb)=z(I);
%         end
%     end
%     end
% end
% 
% 
% % flag low beach bad slope points going north to south
% 
% for xd=min(x):max(x)
%     ndx=find(x == xd);
%     for nys=ndx
%     for yd=max(y(nys)):-1:max(y(nys))-BaseDepth+1
%         mdy1=find( y == yd & x == xd );
%         mdy2=find( y == yd - 1 & x == xd );
%         if z(mdy1) < 1 & ( isempty(mdy2) | z(mdy1) > z(mdy2) )  
%             I=mdy1;
%             nb=nb+1;
%             QC(ldx(CurrentSurveyNumber)).X(nb)=x(I);
%             QC(ldx(CurrentSurveyNumber)).Y(nb)=y(I);
%             QC(ldx(CurrentSurveyNumber)).Z(nb)=z(I);
%         end
%     end
%     end
% end

if numel(QC(ldx(CurrentSurveyNumber)).X)
% remove any duplicate deleted points
% check for duplicate x, y, z
[uxy,ia,ic]=unique([[QC(ldx(CurrentSurveyNumber)).X]',[QC(ldx(CurrentSurveyNumber)).Y]'],'rows');
% keep unique points
QC(ldx(CurrentSurveyNumber)).X=QC(ldx(CurrentSurveyNumber)).X(ia);
QC(ldx(CurrentSurveyNumber)).Y=QC(ldx(CurrentSurveyNumber)).Y(ia);
QC(ldx(CurrentSurveyNumber)).Z=QC(ldx(CurrentSurveyNumber)).Z(ia);
end
    
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

%SelectAreaSA
% loop around

% delete(hpop)
% delete(htxt)
% DeleteAreaSA
