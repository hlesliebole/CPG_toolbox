% Code to create a table of Mop reference frame beach widths from 
% CoastSat (Vos) transect shoreline positions.

% Uses CoastSat transect definitions that have been parsed from 
% the CoastSat_transect_layer.geojson file and placed in the struct
% array "VosT" stored in VosTransects.mat
%
%  VosT also identifies the Mop nearest a CoastSat transect and the
%  xshore correction from the shoreline distance to be relative to the Mop 
%  back beach line.

clearvars

% set the Mop range to include in the table
Mop1=497;
Mop2=682;

% create table with 6 columns of Data 
DateNum=[]; % image datenum
Year=[]; % image year
Month=[]; % image month
Day=[]; % image day
MopNum=[]; % Mop number
BackMopXoffset=[]; % Vos to Mop xshore ref frame offset BackMopDist=VosDist+BackMopXoffset
BackMopDist=[]; % Mop xshore distance (BackMopXoffset correction applied)

% initialize the Mopified Vos shoreline table elements
VosShoreline=table(MopNum,BackMopXoffset,DateNum,Year,Month,Day,BackMopDist);

load VosTransects.mat % load CoastSat transect info
for n=1:size(VosT,2)
    if isempty(VosT(n).BackMop) % fill any empty fields with NaNs
       VosT(n).BackMop=NaN;
       VosT(n).BackMopXoffset=NaN;
    end
end

% loop through Mops, find nearest CoastSat transect and convert
%  shoreline distance values to Mop back beach distances.
for MopNumber=Mop1:Mop2
    idx=find(round([VosT.BackMop]) == MopNumber); 
    if numel(idx) > 0
        idx=idx(1);
    fprintf('%i  %s\n',MopNumber,VosT(idx).Name)
 % fetch CoastSat shoreline data data
  url=['http://coastsat.wrl.unsw.edu.au/time-series/' VosT(idx).Name '/'];
  s=webread(url);

ss=strsplit(s,'\n');
k=0;
for n=1:size(ss,2)
    if ~isempty(ss{n})
       st=strsplit(regexprep(ss{n},',',' '),' ');
       k=k+1;
       
     % add row to shoreline table based on mops with Mhw contour location
     Shoreline.MopNum=MopNumber;
     Shoreline.BackMopXoffset=VosT(idx).BackMopXoffset;
     Shoreline.DateNum=datenum(st{1});
     Shoreline.Year=year(datetime(st{1}));
     Shoreline.Month=month(datetime(st{1}));
     Shoreline.Day=day(datetime(st{1})); 
     Shoreline.BackMopDist=str2double(st{3})+VosT(idx).BackMopXoffset;
     VosShoreline = [VosShoreline;struct2table(Shoreline)];
     
    end   
end
else
    fprintf('No Vos Transect is closest to Mop: %i\n',MopNumber)
end

end

save VosShoreline.mat VosShoreline
