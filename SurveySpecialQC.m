function [SA]=SurveySpecialQC(SA)
%-----------------------------------------------------------------
% special QC section where specific undesired fixed features can be
%  flagged from the survey data set for sediment analysis purposes. 
%  "bad" points have their Class parameter set = -999

%  Add additional code here, using the Cardiff wastewater pipe 
%  as an example, if other features need removal.
%-----------------------------------------------------------------

% loop through surveys
for n=1:size(SA,2)

% Cardiff State Beach Special QC
    % check if in Cardiff area
    if SA(n).Mopnum > 663 && SA(n).Mopnum < 683
        
% convert utm to lat lon
[y,x]=utm2deg(SA(n).X,SA(n).Y,repmat('11 S',[length(SA(n).X) 1]));

%  Cardiff #1 Wastewater Pipe Survey Data Removal

% Remove any survey points that cross over the
% exposed Cardiff wastewater ocean outfall pipeline in > 7m water depth 
ipipe=find(x < -117.2828 & y > 33.0091 & y < 33.0095);
if(~isempty(ipipe)) % if points found in shore box data, remove them
    fprintf('%s Flagging %g Exposed Ocean Outfall Pipe Survey Points.\n',...
        datestr(SA(n).Datenum),length(ipipe))
     SA(n).Class(ipipe)=-999;  
%x(ipipe)=[];y(ipipe)=[];z(ipipe)=[];t(ipipe)=[];
end

% Cardiff #2 has couple bad data points in the 20-Apr-2016 survey 
%  that put a pit in the lowest elevation surface
if SA(n).Datenum == datenum(2016,4,20)

ipit=find(x > -117.2808 & x < -117.280 & y > 33.0083 & y < 33.0091);
if(~isempty(ipit))
fprintf(...
'%s Flagging %g 4/20/16 -z outliers that impact the lowest observed surface.\n',...
    datestr(SA(n).Datenum),length(ipit))
%x(ipit)=[];y(ipit)=[];z(ipit)=[];t(ipit)=[];
SA(n).Class(ipit)=-999; 
end

end

% Cardiff #3: some bad data points in the 18-Jan-2011 survey 
%  that put a pit in the lowest elevation surface on cardiff reef
if SA(n).Datenum == datenum(2011,1,18)

ipit=find( x > -117.2822 & x < -117.2811 & y > 33.0001 & y < 33.0015);
if(~isempty(ipit))
fprintf(...
'%s Flagging %g 1/18/11 -z outliers that impact the lowest observed surface.\n',...
    datestr(SA(n).Datenum),length(ipit))
%x(ipit)=[];y(ipit)=[];z(ipit)=[];t(ipit)=[];
SA(n).Class(ipit)=-999; 
end

end

end % end cardiff area 

% TP State Beach Special QC
    % Mops 534 and 535 in 7 Feb 2005 jumbo have too deep jetski data
    if SA(n).Mopnum > 533 && SA(n).Mopnum < 536
      if SA(n).Datenum == datenum(2005,2,7)
          ibad=find(SA(n).Z < -1.4);
          SA(n).Class(ibad)=-999;
      end
    end
    
    % Mops 555 and 556 in 4 Apr 2016 jumbo have too deep jetski data
    if SA(n).Mopnum > 554 && SA(n).Mopnum < 557
      if SA(n).Datenum == datenum(2016,4,4)
          ibad=find(SA(n).Z < -1.7);
          SA(n).Class(ibad)=-999;
      end
    end

% TP State Beach Special QC
    % check if in TP area
    if SA(n).Mopnum > 534 && SA(n).Mopnum < 608
        
 %  Torrey bad low elv survey transect in 4 Apr 2016 survey       
    if SA(n).Datenum == datenum(2016,4,4)  
        
% convert utm to lat lon
[y,x]=utm2deg(SA(n).X,SA(n).Y,repmat('11 S',[length(SA(n).X) 1]));
    
ipit=find(y > 32.9033 & y < 32.9052);
if(~isempty(ipit))
fprintf(...
'Flagging %g 4/4/16 -z outliers that impact the lowest observed surface.\n',...
    length(ipit))
%x(ipit)=[];y(ipit)=[];z(ipit)=[];t(ipit)=[];
SA(n).Class(ipit)=-999; 
end

    end
    
    end  
    
end    
end