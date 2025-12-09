function [Zmedian]=GetSAmedianProfile(SA)

% returns a global median profile by calculating 
%  monthly->quarterly->seasonal->global medians  
%  of the entire data set.

dtime=datetime([SA.Datenum],'convertfrom','datenum');
y=year(dtime);
mn=month(dtime);
dy=day(dtime);

% month medians
MonMedians(1:12)=NaN;
for m=1:12
  idx=find(mn == m);
%   if ~isempty(idx)
%       MonMedians(m)=median(SA(idx).Z);
      
end