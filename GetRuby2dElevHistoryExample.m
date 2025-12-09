%% example code and function to get all the CPG MOPs survey elevations
%  for a single location during ruby2d.

addpath /volumes/group/MOPS
addpath /volumes/group/MOPS/toolbox

% example TP lat, lon input but passing function x,y in utm is ok too.
y=32.928200;x=-117.260000; 

[SurvDatenum,SurvSrc,NearestPointDist,Z]=GetRuby2dElevHistory(x,y);

% print an example table of returned info
fprintf('\n                      Gridded     Nearest Survey Point\n')
fprintf('SurveyDate    Type   z(m,navd88)      Distance (m)\n\n')

for i=1:numel(SurvDatenum)
    fprintf('%s %8s   %6.2f           %6.1f\n',...
        datestr(SurvDatenum(i)),SurvSrc{i},Z(i),NearestPointDist(i)) 
end

%% ----------------------------------------------------------------------
function [SurvDatenum,SurvSrc,NearestPointDist,Z]=GetRuby2dElevHistory(x,y)

%  returns the survey elevation history of a single point
%  during the Ruby2d experiment Oct 2021-Feb 2022.

%  x,y can either be Eutm,Nutm or Long,Lat

%    SurvDatenum= vector of the survey dates
%    SurvSrc=cell array of the survey types
%    NearestPointDist=vector of distances in meters to the nearest true
%                     survey point
%

% figure out if it is utm or long,lat
%    if long,lat convert to utm
if x < 0
    [x,y,utmzone]=deg2utm(y,x);
end

% find nearest mop area of the points
load MopTableUTM.mat
Nmop=XY2MopNumsV2(x,y,Mop);
%fprintf('Nearest Mop: %i\n', Nmop)

% load SA and SG mat file for this mop
eval(['load ' ['M' num2str(Nmop,'%5.5i') 'SA.mat'] ' SA']);
eval(['load ' ['M' num2str(Nmop,'%5.5i') 'SG.mat'] ' SG']);

% find SA and SG gridded surveys during ruby2d 
%   Skip any RTKdrone data for now.

rgi=find([SG.Datenum] >= datenum(2021,10,1) &...
    [SG.Datenum] < datenum(2022,3,1) & strcmp({SG.Source},'RTKdrone') == 0);


Z=[];
SurvDatenum=[];
SurvSrc={};
NearestPointDist=[];
% loop through gridded surveys
for n=rgi
   
    % find 1m gridded point
    ipt=find( SG(n).X == round(x) & SG(n).Y == round(y) );
    
    % skip survey if no exact gridded point found
    if ~isempty(ipt)
        SurvDatenum=[SurvDatenum SG(n).Datenum];
        SurvSrc=[SurvSrc,{SG(n).Source}];
        Z=[Z SG(n).Z(ipt)];
        % find matching original 1m average Survey data set
        isrv=find([SA.Datenum] == SG(n).Datenum & strcmp({SA.Source},SG(n).Source) == 1);
        % find the nearest pregridded point to the input x,y location
        [k,dist]=dsearchn(horzcat(SA(isrv).X,SA(isrv).Y),[x y]);
        NearestPointDist=[NearestPointDist dist];
    end
end          
     
end
