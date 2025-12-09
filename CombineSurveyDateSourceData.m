function [cx,cy,cz,cclass]=CombineSurveyDateSourceData(SA,SurvNum)

%  
%  Combines airborne lidar survey data for SA struct entries that are 
%  on the same date and from the same source.
%
%  SA(SurvNum) is the SA entry that you want to combine with any
%   other same airborne source survey file entries in SA.
%
%  Returns combined x,y,z, class data.  If no additional surveys are found
%  then these are just the same points as SA(SurvNum).X , .Y , .Z and
%  .Class
%
%  Code will combine same-source data from the following airborne 
%  lidar sources:
%                    CCC,KMair,USACE,USGS,UTAir

cx=SA(SurvNum).X;
cy=SA(SurvNum).Y;
cz=SA(SurvNum).Z;
cclass=SA(SurvNum).Class;

if strcmp(SA(SurvNum).Source,'CCC') == 1 || ... 
        strcmp(SA(SurvNum).Source,'KMair') == 1 || ... 
        strcmp(SA(SurvNum).Source,'USACE') == 1 || ... 
        strcmp(SA(SurvNum).Source,'USGS') == 1 || ... 
        strcmp(SA(SurvNum).Source,'UTAir') == 1 
          
%fprintf('Looking for Same Date/Source Data: %s\n for: %s\n',...
%    SA(SurvNum).Source,SA(SurvNum).File);
        
% find any other entries that have both the same
%  date and the same survey source
      idx=find( [SA.Datenum] == SA(SurvNum).Datenum & ...
          strcmpi({SA.Source}, SA(SurvNum).Source)==1 );
      idx(idx == SurvNum)=[]; % exclude input survey number

% if duplicate date/source entries found, combine
%  their x,y,z data into a single entry and add
%  to list of duplicate entries to delete
      if numel(idx) > 0
          %fprintf('Duplicate Source: %s ; Indices: %s\n',SA(n).Source,num2str(idx))
          
          for i=1:numel(idx)
              fprintf('Combining %i pts\nfrom: %s\n',numel(SA(idx(i)).X),SA(idx(i)).File)
              cx=[cx' SA(idx(i)).X']';
              cy=[cy' SA(idx(i)).Y']';
              cz=[cz' SA(idx(i)).Z']';
              cclass=[cclass' SA(idx(i)).Class']';
          end
          
          % remove any duplicate points
          % check for duplicate x, y with different z
          nxyzpts=size(cx,1);
          [uxy,ia,~]=unique([cx, cy],'rows');
          uniqxy=size(uxy,1);dupxy=nxyzpts-uniqxy;
          fprintf('Removing %i duplicate combined points.\n',dupxy);
          % remove those points 
          cx=cx(ia);cy=cy(ia);cz=cz(ia);cclass=cclass(ia);
          
      end

else
    %fprintf('Code will not combine data from source: %s\n',SA(SurvNum).Source);
end

end





 