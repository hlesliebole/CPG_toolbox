function [BeachWidthMin,BeachWidthMax]=GetCpgProfileBeachWidths(Znavd88,X1D,Z1D) 

% calls the function intersections.m

% Calculates the min and max cross-shore X1D locatiosn that intersect Znavd88
% If there is a sinlge intersection BeachWidthMin = BeachWidthMax.  If
% there is no intersection they are NaNs

% loop through profiles and find Znavd88 intersection with X1D as the
% definition of beach width

nsurv=size(Z1D,1);
BeachWidthMin=NaN(1,nsurv);
BeachWidthMax=NaN(1,nsurv);

for n=1:nsurv
 igood=find(~isnan(Z1D(n,:))); % profile points with valida data

 if numel(igood) > 1 % need at least 2 good profile points

 % find x intersections with Znavd88
 xBW=intersections([X1D(1) X1D(end)],[Znavd88 Znavd88],X1D(igood),Z1D(n,igood));
 
   if ~isempty(xBW)
     % keep the nearest and furthest intersection from the back beach
     % (often the same but useful for catching weird profiles
     BeachWidthMin(n)=xBW(1); 
     BeachWidthMax(n)=xBW(end);
   else
     BeachWidthMin(n)=NaN;
     BeachWidthMax(n)=NaN;
   end

 else
     BeachWidthMin(n)=NaN;
     BeachWidthMax(n)=NaN;
 end

end

end
