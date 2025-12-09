%% Example code to plot a composite of 2009 and 2014 SHOALS 1m spatial averaged survey 
%   data for a reach of Mop points. Shoals surveys occurred over mutiple days in 2009 and 2014.

%% Combine 1m avg Mop SA struct arrays for a specified (starting,ending) Mop reach 

%CS=SAcombineMops(481,506);  % la jolla cove to the marine room
%CS=SAcombineMops(665,715);  % cardiff SB to north of encinitas point/swamis
%CS=SAcombineMops(658,686);  % cardiff SB 
CS=SAcombineMops(615,640);  % Del Mar

%% find the USACE shoals surveys in the combined CS survey struct array

%% option for just the 2009 or 2014 survey dates
% iyr=2009;
% iyr=2014;
% sidx=find(strcmp({CS.Source},'USACE') & year(datetime([CS.Datenum],'convertfrom','datenum')) == iyr);%%  

%% option to find all of the shoals surveys in both years 
sidx=find(strcmp({CS.Source},'USACE'));

%% get global x,y limits for the composite grid
xmin=[];xmax=[];ymin=[];ymax=[];
for n=sidx  % loop through shoals surveys
 xmin=min([xmin [CS(n).X]']);xmax=max([xmax [CS(n).X]']);
 ymin=min([ymin [CS(n).Y]']);ymax=max([ymax [CS(n).Y]']);
end

%% make 2d X,Y grid arrays
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);
% initialize composite elev grid as no data NaNs
Z=nan(size(X)); 

%% progressively overlay Shoals 1m average data onto the grid 
%  so any grid points with redundant info keep the most recent survey
for n=sidx % loop through shoal survey dates
    % 1d grid indices with data
    gdx=sub2ind(size(X),round(CS(n).Y)-ymin+1,round(CS(n).X)-xmin+1); 
    Z(gdx)=CS(n).Z; % overlay on grid
end

%% reduce to x,y,z vectors of just the composite grid points with data
ndx=find(~isnan(Z(:)));
x=X(ndx);y=Y(ndx);z=Z(ndx);

%% Make some plots

% 2D plot

figure('position',[143   264   811   523]);
%ScatterPlotBeachUTM(x,y,z,'2d');BeachColorbar; % beach color option
ScatterPlotSubaqueousUTM(x,y,z,-12,'2d'); % subaqueous only with max colors
title([{'USACE SHOALS Airborne Topo-Bathymetric Surveys:'};...
    {[ datestr(CS(sidx(1)).Datenum) ' to ' datestr(CS(sidx(end)).Datenum)]}])

PlotLabelMopTransectUTM(667,'2d','m','ShadeOff')
PlotLabelMopTransectUTM(669,'2d','m','ShadeOff')

% 3D plot

figure('position',[243   164   811   523]);
%ScatterPlotBeachUTM(x,y,z,'3d');BeachColorbar; % beach color option
ScatterPlotSubaqueousUTM(x,y,z,-12,'3d'); % subaqueous only with max colors
title([{'USACE SHOALS Airborne Topo-Bathymetric Surveys:'};...
    {[ datestr(CS(sidx(1)).Datenum) ' to ' datestr(CS(sidx(end)).Datenum)]}])

PlotLabelMopTransectUTM(667,'2d','m','ShadeOff')
PlotLabelMopTransectUTM(669,'2d','m','ShadeOff')
view(-127.5,30)
