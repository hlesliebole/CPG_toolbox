function [BeachVolMin,BeachVolMax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1D,Z1D)


% calls the function intersections.m

% Calculates the min and max cross-shore X1D locatiosn that intersect Znavd88
% If there is a sinlge intersection BeachWidthMin = BeachWidthMax.  If
% there is no intersection they are NaNs

% loop through profiles and find Znavd88 intersection with X1D as the
% definition of beach width

nsurv=size(Z1D,1);
BeachVolMin=NaN(1,nsurv);
BeachVolMax=NaN(1,nsurv);
BackGapSize=NaN(1,nsurv);

for n=1:nsurv

 % find data gap sizes
 % identify the size of the gap each xshore point is in (0 = not in a gap)
     sz=gapsize(Z1D(n,:));
     idxlast=find(sz == 0,1,'first'); % first xshore valid data point index
     idxback=find(X1D == 0); % back beach x=0 point index
     if ~isempty(idxlast-idxback)
      BackGapSize(n)=idxlast-idxback; % gap between back beach point and first data point 
     else
      BackGapSize(n)=NaN;
     end

     % set gap points equal to first valid data point value if gap <=
     % XbackgapTol
     %fprintf('%i back gap %i\n',n,BackGapSize)
     if BackGapSize(n) > 0 & BackGapSize(n) <= XbackgapTol
         for idx=idxback:idxlast-1
             Z1D(n,idx)=Z1D(n,idxlast);
             sz(idx)=0; % set gap size to 0 for this point
         end
     end


 igood=find(~isnan(Z1D(n,:))); % profile points with valida data

 if numel(igood) > 1 % need at least 2 good profile points
 % find x intersections with input Znavd88 elevation
 xBW=intersections([X1D(1) X1D(end)],[Znavd88 Znavd88],X1D(igood),Z1D(n,igood));
 %fprintf('%i intersect %f8.1 %f8.1\n',n,xBW(1),xBW(end))
 

 if ~isempty(xBW)
     xBW=round(xBW); % round intersections to 1m xshore resolution indices

     % set any profile points below the volume elevation cutoff to = the
     %  cutoff (ie. no negative volumes when summing volume above the elev)
     ilow=find(Z1D(n,:) < Znavd88);
     if ~isempty(ilow)
         Z1D(n,ilow)=Znavd88;
     end

 % find if there are any gaps between the nearest xshore intersect and the back
 % beach
 idxBW=find(X1D == xBW(1));

 if idxBW > idxback
     %fprintf('Min %i %i\n',n,sum(sz(idxback:idxBW)) )
    if sum( sz(idxback:idxBW)) == 0
 % if no gaps, calculate a volume
     BeachVolMin(n)=sum(Z1D(n,idxback:idxBW)-Znavd88); 
    end
   end

 % find if there are any gaps between the farthest xshore intersect and the back
 % beach
 idxBW=find(X1D == xBW(end));

 if idxBW > idxback
     %fprintf('Min %i %i\n',n,sum(sz(idxback:idxBW)) )
    if sum( sz(idxback:idxBW)) == 0
 % if no gaps, calculate a volume
     BeachVolMax(n)=sum(Z1D(n,idxback:idxBW)-Znavd88); 
    end
   end

 end % end if intersection was found
 end % end if enough data to look for an intersection

end

 fprintf('Median back gap size: %i\n',median(BackGapSize,'omitnan'))
 fprintf('Number of valid Volumes: %i/%i\n',...
     numel(find(~isnan(BeachVolMax))),nsurv)




end
