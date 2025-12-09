%function [xgimg,ygimg,sbgimg]=getShoreBoxGoogleMap(Mop,mop1,mop2)

% gets google map and transforms into shore box coordinates
%
% eg. for Cardiff
%
%   [xgimg,ygimg,sbgimg]=getShoreBoxGoogleMap(Mop,666,683);
%   image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');
% 
% input
%
%  Mop: mop info table (load from MopTable.mat)
%
%  mop1: starting numeric mop number (all-CA numbering convention)
%  mop2: ending numeric mop number (all-CA numbering convention)
%
% returns:
%
%  xgimag: shore box x coords of image
%  ygimag: shore box y coords of image
%  sbgimag: google map image transformed ito shore box coords

%
%
%      
%
%
%  The map area is edgefactor*100%% larger than than the shore box on 
%  all sides to allow for displaying more of the landward image if desired.
addpath /Users/William/Desktop/Ron

load MopTable

% mop1=1140;
% mop2=1210;
mop1=860;
mop2=916;
mop1=1211;
mop2=1316;
mop1=755; %Carlsbad
mop2=870;
mop1=512; % Blacks to Oside harbor
mop2=929;
mop1=594; % Del Mar
mop2=646;
mop1=576; % North Torrey
mop2=598;
mop1=635; % Solana Cardiff
mop2=685;
mop1=493; % LJ Shores
mop2=517;
mop1=635; % Solana 
mop2=665;


edgefactor=0.1;

% get google map for map area that contains the specified mop
%  range lines
x1=min([Mop.BackLon(mop1:mop2)' Mop.OffLon(mop1:mop2)']); 
y1=min([Mop.BackLat(mop1:mop2)' Mop.OffLat(mop1:mop2)']);
x2=max([Mop.BackLon(mop1:mop2)' Mop.OffLon(mop1:mop2)']); 
y2=max([Mop.BackLat(mop1:mop2)' Mop.OffLat(mop1:mop2)']); 

gf=figure('position',[150    61   844   733]);
gax=axes('xlim',[x1 x2],'ylim',[y1 y2]);
plot_google_map('MapType', 'satellite','Alpha', 1,'axis',gax);
% option that just returns plotting info
[lonv,latv,img]=plot_google_map('MapType', 'satellite','Alpha', 1,'axis',gax);
%close(gf);

fprintf('Transforming the Google Map.  Takes a few minutes...\n')
[X2d,Y2d] = meshgrid(lonv,latv);
% transform google image coords to shore box coords
[MopSB,xsb,ysb]=llz2shorebox(Mop(mop1:mop2,:),X2d(:),Y2d(:));

% get shore box coord range based on transformed mop transect range
xsb1=round(min([MopSB.BackLon' MopSB.OffLon'])); 
ysb1=round(min([MopSB.BackLat' MopSB.OffLat']));
xsb2=round(max([MopSB.BackLon' MopSB.OffLon'])); 
ysb2=round(max([MopSB.BackLat' MopSB.OffLat'])); 
% clip out shore box area + (edgefactor*100)% on all sides (desired map view area)
dx=round(edgefactor*(xsb2-xsb1));
dy=round(edgefactor*(ysb2-ysb1));
xi=xsb1-dx:xsb2+dx;
yi=ysb1-dy*5:ysb2+dy*3;
% find transformed google image indices in shorebox view area
i=find(round(xsb) >= xi(1) & round(xsb) <= xi(end) & ...
    round(ysb) >= yi(1) & round(ysb) <= yi(end));

r=img(:,:,1);r=r(i);
b=img(:,:,2);b=b(i);
g=img(:,:,3);g=g(i);

% create empty shorebox image matrix
sbimg=NaN(length(yi),length(xi),3);
% grid transfomed image r g b values into regular grid
sbimg(:,:,1)=griddata(double(xsb(i)),double(ysb(i)'),double(r),double(xi),double(yi'));
sbimg(:,:,2)=griddata(double(xsb(i)),double(ysb(i)'),double(b),double(xi),double(yi'));
sbimg(:,:,3)=griddata(double(xsb(i)),double(ysb(i)'),double(g),double(xi),double(yi'));
% define image variables to be returned
sbgimg=uint8(sbimg);
xgimg=xi;
ygimg=yi';

%image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');
%save SanOnofreShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save OceansideShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save SanClementeShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save CarlsbadShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save OsideCellShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save DelMarShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save TorreyShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save SolanaCardiffShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
%save LJShoresShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg
save SolanaShoreboxMap.mat MopSB xsb ysb xgimg ygimg sbgimg


%end

%----------------------Plot Google Map----------------------------

% Coordinate transformation functions

function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);
end

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;
end

function ZI = myTurboInterp2(X,Y,Z,XI,YI)
% An extremely fast nearest neighbour 2D interpolation, assuming both input
% and output grids consist only of squares, meaning:
% - uniform X for each column
% - uniform Y for each row
XI = XI(1,:);
X = X(1,:);
YI = YI(:,1);
Y = Y(:,1);

xiPos = nan*ones(size(XI));
xLen = length(X);
yiPos = nan*ones(size(YI));
yLen = length(Y);
% find x conversion
xPos = 1;
for idx = 1:length(xiPos)
    if XI(idx) >= X(1) && XI(idx) <= X(end)
        while xPos < xLen && X(xPos+1)<XI(idx)
            xPos = xPos + 1;
        end
        diffs = abs(X(xPos:xPos+1)-XI(idx));
        if diffs(1) < diffs(2)
            xiPos(idx) = xPos;
        else
            xiPos(idx) = xPos + 1;
        end
    end
end
% find y conversion
yPos = 1;
for idx = 1:length(yiPos)
    if YI(idx) <= Y(1) && YI(idx) >= Y(end)
        while yPos < yLen && Y(yPos+1)>YI(idx)
            yPos = yPos + 1;
        end
        diffs = abs(Y(yPos:yPos+1)-YI(idx));
        if diffs(1) < diffs(2)
            yiPos(idx) = yPos;
        else
            yiPos(idx) = yPos + 1;
        end
    end
end
ZI = Z(yiPos,xiPos,:);
end

function update_google_map(obj,evd)
% callback function for auto-refresh
drawnow;
try
    axHandle = evd.Axes;
catch ex
    % Event doesn't contain the correct axes. Panic!
    axHandle = gca;
end
ud = get(axHandle, 'UserData');
if isfield(ud, 'gmap_params')
    params = ud.gmap_params;
    plot_google_map(params{:});
end
end

function update_google_map_fig(obj,evd)
% callback function for auto-refresh
drawnow;
axes_objs = findobj(get(gcf,'children'),'type','axes');
for idx = 1:length(axes_objs)
    if ~isempty(findobj(get(axes_objs(idx),'children'),'tag','gmap'));
        ud = get(axes_objs(idx), 'UserData');
        if isfield(ud, 'gmap_params')
            params = ud.gmap_params;
        else
            params = {};
        end
        
        % Add axes to inputs if needed
        if ~sum(strcmpi(params, 'Axis'))
            params = [params, {'Axis', axes_objs(idx)}];
        end
        plot_google_map(params{:});
    end
end
end

function cleanupFunc(h)
ud = get(h, 'UserData');
if isstruct(ud) && isfield(ud, 'gmap_params')
    ud = rmfield(ud, 'gmap_params');
    set(h, 'UserData', ud);
end
end

function varargout = plot_google_map(varargin)
% function h = plot_google_map(varargin)
% Plots a google map on the current axes using the Google Static Maps API
%
% USAGE:
% h = plot_google_map(Property, Value,...)
% Plots the map on the given axes. Used also if no output is specified
%
% Or:
% [lonVec latVec imag] = plot_google_map(Property, Value,...)
% Returns the map without plotting it
%
% PROPERTIES:
%    Axis           - Axis handle. If not given, gca is used.
%    Height (640)   - Height of the image in pixels (max 640)
%    Width  (640)   - Width of the image in pixels (max 640)
%    Scale (2)      - (1/2) Resolution scale factor. Using Scale=2 will
%                     double the resulotion of the downloaded image (up
%                     to 1280x1280) and will result in finer rendering,
%                     but processing time will be longer.
%    Resize (1)     - (recommended 1-2) Resolution upsampling factor. 
%                     Increases image resolution using imresize(). This results
%                     in a finer image but it needs the image processing
%                     toolbox and processing time will be longer.
%    MapType        - ('roadmap') Type of map to return. Any of [roadmap, 
%                     satellite, terrain, hybrid]. See the Google Maps API for
%                     more information. 
%    Alpha (1)      - (0-1) Transparency level of the map (0 is fully
%                     transparent). While the map is always moved to the
%                     bottom of the plot (i.e. will not hide previously
%                     drawn items), this can be useful in order to increase
%                     readability if many colors are plotted 
%                     (using SCATTER for example).
%    ShowLabels (1) - (0/1) Controls whether to display city/street textual labels on the map
%    Style          - (string) A style configuration string. See:
%                     https://developers.google.com/maps/documentation/static-maps/?csw=1#StyledMaps
%                     http://instrument.github.io/styled-maps-wizard/
%    Language       - (string) A 2 letter ISO 639-1 language code for displaying labels in a 
%                     local language instead of English (where available).
%                     For example, for Chinese use:
%                     plot_google_map('language','zh')
%                     For the list of codes, see:
%                     http://en.wikipedia.org/wiki/List_of_ISO_639-1_codes
%    Marker         - The marker argument is a text string with fields
%                     conforming to the Google Maps API. The
%                     following are valid examples:
%                     '43.0738740,-70.713993' (default midsize orange marker)
%                     '43.0738740,-70.713993,blue' (midsize blue marker)
%                     '43.0738740,-70.713993,yellowa' (midsize yellow
%                     marker with label "A")
%                     '43.0738740,-70.713993,tinyredb' (tiny red marker
%                     with label "B")
%    Refresh (1)    - (0/1) defines whether to automatically refresh the
%                     map upon zoom/pan action on the figure.
%    AutoAxis (1)   - (0/1) defines whether to automatically adjust the axis
%                     of the plot to avoid the map being stretched.
%                     This will adjust the span to be correct
%                     according to the shape of the map axes.
%    MapScale (0)  - (0/1) defines wheteher to add a scale indicator to
%                     the map.
%    ScaleWidth (0.25) - (0.1-0.9) defines the max width of the scale
%                     indicator relative to the map width.
%    ScaleLocation (sw) - (ne, n, se, s, sw, nw) defines the location of
%                     scale indicator on the map.
%    ScaleUnits (si) - (si/imp) changes the scale indicator units between 
%                     SI and imperial units.
%    FigureResizeUpdate (1) - (0/1) defines whether to automatically refresh the
%                     map upon resizing the figure. This will ensure map
%                     isn't stretched after figure resize.
%    APIKey         - (string) set your own API key which you obtained from Google: 
%                     http://developers.google.com/maps/documentation/staticmaps/#api_key
%                     This will enable up to 25,000 map requests per day, 
%                     compared to a few hundred requests without a key. 
%                     To set the key, use:
%                     plot_google_map('APIKey','SomeLongStringObtaindFromGoogle')
%                     You need to do this only once to set the key.
%                     To disable the use of a key, use:
%                     plot_google_map('APIKey','')
%
% OUTPUT:
%    h              - Handle to the plotted map
%
%    lonVect        - Vector of Longidute coordinates (WGS84) of the image 
%    latVect        - Vector of Latidute coordinates (WGS84) of the image 
%    imag           - Image matrix (height,width,3) of the map
%
% EXAMPLE - plot a map showing some capitals in Europe:
%    lat = [48.8708   51.5188   41.9260   40.4312   52.523   37.982];
%    lon = [2.4131    -0.1300    12.4951   -3.6788    13.415   23.715];
%    plot(lon, lat, '.r', 'MarkerSize', 20)
%    plot_google_map('MapScale', 1)
%
% References:
%  http://www.mathworks.com/matlabcentral/fileexchange/24113
%  http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
%  http://developers.google.com/maps/documentation/staticmaps/
%  https://www.mathworks.com/matlabcentral/fileexchange/33545-automatic-map-scale-generation
%
% Acknowledgements:
%  Val Schmidt for the submission of get_google_map.m
%  Jonathan Sullivan for the submission of makescale.m
%
% Author:
%  Zohar Bar-Yehuda
%
% Version 2.0 - 08/04/2018
%       - Add an option to show a map scale
%       - Several bugfixes
% Version 1.8 - 25/04/2016 - By Hannes Diethelm
%       - Add resize parameter to resize image using imresize()
%       - Fix scale parameter
% Version 1.7 - 14/04/2016
%       - Add custom style support
% Version 1.6 - 12/11/2015
%       - Use system temp folder for writing image files (with fallback to current dir if missing write permissions)
% Version 1.5 - 20/11/2014
%       - Support for MATLAB R2014b
%       - several fixes for complex layouts: several maps in one figure, 
%         map inside a panel, specifying axis handle as input (thanks to Luke Plausin)
% Version 1.4 - 25/03/2014
%       - Added the language parameter for showing labels in a local language
%       - Display the URL on error to allow easier debugging of API errors
% Version 1.3 - 06/10/2013
%       - Improved functionality of AutoAxis, which now handles any shape of map axes. 
%         Now also updates the extent of the map if the figure is resized.
%       - Added the showLabels parameter which allows hiding the textual labels on the map.
% Version 1.2 - 16/06/2012
%       - Support use of the "scale=2" parameter by default for finer rendering (set scale=1 if too slow).
%       - Auto-adjust axis extent so the map isn't stretched.
%       - Set and use an API key which enables a much higher usage volume per day.
% Version 1.1 - 25/08/2011

persistent apiKey useTemp
apiKey = google_api_key_local();

if isempty(useTemp)
    % first run, check if we have wrtie access to the temp folder
    try 
        tempfilename = tempname;
        fid = fopen(tempfilename, 'w');
        if fid > 0
            fclose(fid);
            useTemp = true;
            delete(tempfilename);
        else
            % Don't have write access to temp folder or it doesn't exist, fallback to current dir
            useTemp = false;
        end
    catch
        % in case tempname fails for some reason
        useTemp = false;
    end
end

hold on

% Default parametrs
axHandle = gca;
set(axHandle, 'Layer','top'); % Put axis on top of image, so it doesn't hide the axis lines and ticks
height = 640;
width = 640;
scale = 2;
resize = 1;
maptype = 'roadmap';
alphaData = 1;
autoRefresh = 1;
figureResizeUpdate = 1;
autoAxis = 1;
showLabels = 1;
language = '';
markeridx = 1;
markerlist = {};
style = '';
mapScale = 0;
scaleWidth = 0.25;
scaleLocation = 'se';
scaleUnits = 'si';

% Handle input arguments
if nargin >= 2
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'axis'
                axHandle = varargin{idx+1};
            case 'height'
                height = varargin{idx+1};
            case 'width'
                width = varargin{idx+1};
            case 'scale'
                scale = round(varargin{idx+1});
                if scale < 1 || scale > 2
                    error('Scale must be 1 or 2');
                end
            case 'resize'
                resize = varargin{idx+1};
            case 'maptype'
                maptype = varargin{idx+1};
            case 'alpha'
                alphaData = varargin{idx+1};
            case 'refresh'
                autoRefresh = varargin{idx+1};
            case 'showlabels'
                showLabels = varargin{idx+1};
            case 'figureresizeupdate'
                figureResizeUpdate = varargin{idx+1};
            case 'language'
                language = varargin{idx+1};
            case 'marker'
                markerlist{markeridx} = varargin{idx+1};
                markeridx = markeridx + 1;
            case 'autoaxis'
                autoAxis = varargin{idx+1};
            case 'apikey'
                apiKey = varargin{idx+1}; % set new key
                % save key to file
                funcFile = which('plot_google_map.m');
                pth = fileparts(funcFile);
                keyFile = fullfile(pth,'api_key.mat');
                save(keyFile,'apiKey')
            case 'style'
                style = varargin{idx+1};
            case 'mapscale'
                mapScale = varargin{idx+1};
            case 'scalewidth'
                scaleWidth = varargin{idx+1};
            case 'scalelocation'
                scaleLocation = varargin{idx+1};
            case 'scaleunits'
                scaleUnits = varargin{idx+1};
            otherwise
                error(['Unrecognized variable: ' varargin{idx}])
        end
    end
end
if height > 640
    height = 640;
end
if width > 640
    width = 640;
end

% Store paramters in axis handle (for auto refresh callbacks)
ud = get(axHandle, 'UserData');
if isempty(ud)
    % explicitly set as struct to avoid warnings
    ud = struct;
end
ud.gmap_params = varargin;
set(axHandle, 'UserData', ud);

curAxis = axis(axHandle);
if max(abs(curAxis)) > 500 || curAxis(3) > 90 || curAxis(4) < -90
    warning('Axis limits are not reasonable for WGS1984, ignoring. Please make sure your plotted data in WGS1984 coordinates,')
    return;
end    

% Enforce Latitude constraints of EPSG:900913 
if curAxis(3) < -85
    curAxis(3) = -85;
end
if curAxis(4) > 85
    curAxis(4) = 85;
end
% Enforce longitude constrains
if curAxis(1) < -180
    curAxis(1) = -180;
end
if curAxis(1) > 180
    curAxis(1) = 0;
end
if curAxis(2) > 180
    curAxis(2) = 180;
end
if curAxis(2) < -180
    curAxis(2) = 0;
end

if isequal(curAxis,[0 1 0 1]) % probably an empty figure
    % display world map
    curAxis = [-200 200 -85 85];
    axis(curAxis)
end


if autoAxis
    % adjust current axis limit to avoid strectched maps
    [xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
    xExtent = diff(xExtent); % just the size of the span
    yExtent = diff(yExtent); 
    % get axes aspect ratio
    drawnow
    org_units = get(axHandle,'Units');
    set(axHandle,'Units','Pixels')
    ax_position = get(axHandle,'position');        
    set(axHandle,'Units',org_units)
    aspect_ratio = ax_position(4) / ax_position(3);
    
    if xExtent*aspect_ratio > yExtent        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
       
        % enlarge the Y extent
        spanY = spanY*xExtent*aspect_ratio/yExtent; % new span
        if spanY > 85
            spanX = spanX * 85 / spanY;
            spanY = spanY * 85 / spanY;
        end
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    elseif yExtent > xExtent*aspect_ratio
        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
        % enlarge the X extent
        spanX = spanX*yExtent/(xExtent*aspect_ratio); % new span
        if spanX > 180
            spanY = spanY * 180 / spanX;
            spanX = spanX * 180 / spanX;
        end
        
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    end            
    % Enforce Latitude constraints of EPSG:900913
    if curAxis(3) < -85
        curAxis(3:4) = curAxis(3:4) + (-85 - curAxis(3));
    end
    if curAxis(4) > 85
        curAxis(3:4) = curAxis(3:4) + (85 - curAxis(4));
    end
    axis(axHandle, curAxis); % update axis as quickly as possible, before downloading new image
    drawnow
end

% Delete previous map from plot (if exists)
if nargout <= 1 % only if in plotting mode
    curChildren = get(axHandle,'children');
    map_objs = findobj(curChildren,'tag','gmap');
    bd_callback = [];
    for idx = 1:length(map_objs)
        if ~isempty(get(map_objs(idx),'ButtonDownFcn'))
            % copy callback properties from current map
            bd_callback = get(map_objs(idx),'ButtonDownFcn');
        end
    end
    ud = get(axHandle, 'UserData');
    delete(map_objs);
    delete(findobj(curChildren,'tag','MapScale'));
    % Recover userdata of axis (cleared in cleanup function)
    set(axHandle, 'UserData', ud);
end

% Calculate zoom level for current axis limits
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
minResX = diff(xExtent) / width;
minResY = diff(yExtent) / height;
minRes = max([minResX minResY]);
tileSize = 256;
initialResolution = 2 * pi * 6378137 / tileSize; % 156543.03392804062 for tileSize 256 pixels
zoomlevel = floor(log2(initialResolution/minRes));

% Enforce valid zoom levels
if zoomlevel < 0 
    zoomlevel = 0;
end
if zoomlevel > 19 
    zoomlevel = 19;
end

% Calculate center coordinate in WGS1984
lat = (curAxis(3)+curAxis(4))/2;
lon = (curAxis(1)+curAxis(2))/2;

% Construct query URL
preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location = ['?center=' num2str(lat,10) ',' num2str(lon,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];
if ~isempty(apiKey)
    keyStr = ['&key=' apiKey];
else
    keyStr = '';
end
markers = '&markers=';
for idx = 1:length(markerlist)
    if idx < length(markerlist)
        markers = [markers markerlist{idx} '%7C'];
    else
        markers = [markers markerlist{idx}];
    end
end

if showLabels == 0
    if ~isempty(style)
        style = [style '&style='];
    end
    style = [style 'feature:all|element:labels|visibility:off'];
end

if ~isempty(language)
    languageStr = ['&language=' language];
else
    languageStr = '';
end
    
if ismember(maptype,{'satellite','hybrid'})
    filename = 'tmp.jpg';
    format = '&format=jpg';
    convertNeeded = 0;
else
    filename = 'tmp.png';
    format = '&format=png';
    convertNeeded = 1;
end
sensor = '&sensor=false';

if ~isempty(style)
    styleStr = ['&style=' style];
else
    styleStr = '';
end

url = [preamble location zoomStr sizeStr maptypeStr format markers languageStr sensor keyStr styleStr];

% Get the image
if useTemp
    filepath = fullfile(tempdir, filename);
else
    filepath = filename;
end

try
    urlwrite(url,filepath);
catch % error downloading map
    warning(['Unable to download map form Google Servers.\n' ...
        'Matlab error was: %s\n\n' ...
        'Possible reasons: missing write permissions, no network connection, quota exceeded, or some other error.\n' ...
        'Consider using an API key if quota problems persist.\n\n' ...
        'To debug, try pasting the following URL in your browser, which may result in a more informative error:\n%s'], lasterr, url);
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    return
end

[M, Mcolor] = imread(filepath);
Mcolor = uint8(Mcolor * 255);
%M = cast(M,'double');
delete(filepath); % delete temp file
width = size(M,2);
height = size(M,1);

% We now want to convert the image from a colormap image with an uneven
% mesh grid, into an RGB truecolor image with a uniform grid.
% This would enable displaying it with IMAGE, instead of PCOLOR.
% Advantages are:
% 1) faster rendering
% 2) makes it possible to display together with other colormap annotations (PCOLOR, SCATTER etc.)

% Convert image from colormap type to RGB truecolor (if PNG is used)
if convertNeeded
    imag = zeros(height,width,3, 'uint8');
    for idx = 1:3
        cur_map = Mcolor(:,idx);
        imag(:,:,idx) = reshape(cur_map(M+1),height,width);
    end
else
    imag = M;
end
% Resize if needed
if resize ~= 1
    imag = imresize(imag, resize, 'bilinear');
end

% Calculate a meshgrid of pixel coordinates in EPSG:900913
width = size(imag,2);
height = size(imag,1);
centerPixelY = round(height/2);
centerPixelX = round(width/2);
[centerX,centerY] = latLonToMeters(lat, lon ); % center coordinates in EPSG:900913
curResolution = initialResolution / 2^zoomlevel / scale / resize; % meters/pixel (EPSG:900913)
xVec = centerX + ((1:width)-centerPixelX) * curResolution; % x vector
yVec = centerY + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh,yMesh] = meshgrid(xVec,yVec); % construct meshgrid 

% convert meshgrid to WGS1984
[lonMesh,latMesh] = metersToLatLon(xMesh,yMesh);

% Next, project the data into a uniform WGS1984 grid
uniHeight = round(height*resize);
uniWidth = round(width*resize);
latVect = linspace(latMesh(1,1),latMesh(end,1),uniHeight);
lonVect = linspace(lonMesh(1,1),lonMesh(1,end),uniWidth);
[uniLonMesh,uniLatMesh] = meshgrid(lonVect,latVect);
uniImag = zeros(uniHeight,uniWidth,3);

% Fast Interpolation to uniform grid
uniImag =  myTurboInterp2(lonMesh,latMesh,imag,uniLonMesh,uniLatMesh);

if nargout <= 1 % plot map
    % display image
    hold(axHandle, 'on');
    cax = caxis;
    h = image(lonVect,latVect,uniImag, 'Parent', axHandle);
    caxis(cax); % Preserve caxis that is sometimes changed by the call to image()
    set(axHandle,'YDir','Normal')
    set(h,'tag','gmap')
    set(h,'AlphaData',alphaData)
    
    % add a dummy image to allow pan/zoom out to x2 of the image extent
    h_tmp = image(lonVect([1 end]),latVect([1 end]),zeros(2),'Visible','off', 'Parent', axHandle, 'CDataMapping', 'scaled');
    set(h_tmp,'tag','gmap')
   
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    axis(axHandle, curAxis) % restore original zoom
    if nargout == 1
        varargout{1} = h;
    end
    set(h, 'UserData', onCleanup(@() cleanupFunc(axHandle)));
    
    % if auto-refresh mode - override zoom callback to allow autumatic 
    % refresh of map upon zoom actions.
    figHandle = axHandle;
    while ~strcmpi(get(figHandle, 'Type'), 'figure')
        % Recursively search for parent figure in case axes are in a panel
        figHandle = get(figHandle, 'Parent');
    end
    
    zoomHandle = zoom(axHandle);   
    panHandle = pan(figHandle); % This isn't ideal, doesn't work for contained axis    
    if autoRefresh        
        set(zoomHandle,'ActionPostCallback',@update_google_map);          
        set(panHandle, 'ActionPostCallback', @update_google_map);        
    else % disable zoom override
        set(zoomHandle,'ActionPostCallback',[]);
        set(panHandle, 'ActionPostCallback',[]);
    end
    
    % set callback for figure resize function, to update extents if figure
    % is streched.
    if figureResizeUpdate &&isempty(get(figHandle, 'ResizeFcn'))
        % set only if not already set by someone else
        set(figHandle, 'ResizeFcn', @update_google_map_fig);       
    end    
    
    % set callback properties 
    set(h,'ButtonDownFcn',bd_callback);
    
    if mapScale
       makescale(axHandle, 'set_callbacks', 0, 'units', scaleUnits, ...
                 'location', scaleLocation, 'width', scaleWidth);
    end
else % don't plot, only return map
    varargout{1} = lonVect;
    varargout{2} = latVect;
    varargout{3} = uniImag;
end
end

function [MopSB,xsb,ysb]=llz2shorebox(Mop,lon,lat)

%  Returns survey lon,lat points in Mop shorebox coordinates (m)
%  where the shorebox (alongshore) x-axis is defined as the line 
%  connecting the first and last Mop backbeach points.  The x,y origin
%  is defined relative to this line, with the x=1,y=1 is set to 
%  encompass all the MOP back beach points in the defined reach.
%  

% input 
%  Mop: Mop table for the specific mop reach of the shorebox
%  lon: Longitudes of the surrounding area's survey data
%  lat: Latiitudes of the surrounding area's survey data
%
% output 
%  MopSB: Mop table with transformed back beach and offshore point coords
%  xsb: Longitudes transformed to shorebox coords (m)
%  ysb: Latiitudes transformed to shorebox coords (m)
%
% returns mop table and LLZ data (in UTM x,y,z format) for shore box
%  area defined by the mop alongshore transects in the Mop table
%  variable and the specified area bounds around those mops.

 
MopSB=Mop;
% convert mop transect back beach and offshore points to utm
[MopSB.BackLon,MopSB.BackLat,utmzone]=deg2utm(Mop.BackLat,Mop.BackLon);
[MopSB.OffLon,MopSB.OffLat,utmzone]=deg2utm(Mop.OffLat,Mop.OffLon);

% initial axis and (0,0) origin based on first and last Mop
x0=MopSB.BackLon(1);y0=MopSB.BackLat(1);
xe=MopSB.BackLon(end);ye=MopSB.BackLat(end);
theta0=atan2((ye-y0),(xe-x0)); % angle of axis
% transform all the Mop back beach points
B=sqrt((MopSB.BackLon-x0).^2+(MopSB.BackLat-y0).^2);
theta=atan2((MopSB.BackLat-y0),(MopSB.BackLon-x0));
thetat=(theta-theta0); % rotation
MopSB.BackLon=B.*cos(thetat); % rotated mop backbeach x values 
MopSB.BackLat=B.*sin(thetat); % rotated mop backbeach y values
% transform all the Mop offshore points
B=sqrt((MopSB.OffLon-x0).^2+(MopSB.OffLat-y0).^2);
theta=atan2((MopSB.OffLat-y0),(MopSB.OffLon-x0));
thetat=(theta-theta0); % rotation
MopSB.OffLon=B.*cos(thetat); % rotated mop backbeach x values 
MopSB.OffLat=B.*sin(thetat); % rotated mop backbeach y values

% Now make final axis origin by shifting y axis values so mininum mop
% back beach y value = 0

dy=min(MopSB.BackLat);

% apply y shift to MopSB y coords
MopSB.BackLat=MopSB.BackLat-dy;
MopSB.OffLat=MopSB.OffLat-dy;

% now transform survey data

% convert lat lon data to utm xdu ydu
[xdu,ydu,utmzone]=deg2utm(lat,lon);

B=sqrt((xdu-x0).^2+(ydu-y0).^2);
theta=atan2((ydu-y0),(xdu-x0));
thetat=(theta-theta0); % rotation
xsb=B.*cos(thetat); % rotated mop backbeach x values 
ysb=B.*sin(thetat); % rotated mop backbeach y values
% apply y shift to survey y coords
ysb=ysb-dy;

end

% -------------------------------------------------------------------------
function  [x,y,utmzone] = deg2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
% Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.
%
% Inputs:
%    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
%    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84
%
% Outputs:
%    x, y , utmzone.   See example
%
% Example 1:
%    Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
%    Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
%    [x,y,utmzone] = deg2utm(Lat,Lon);
%    fprintf('%7.0f ',x)
%       458731  407653  239027  230253  343898  362850
%    fprintf('%7.0f ',y)
%      4462881 5126290 4163083 3171843 4302285 2772478
%    utmzone =
%       30 T
%       32 T
%       11 S
%       28 R
%       15 S
%       51 R
%
% Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
%    LatDMS=[40 18 55.56; 46 17 2.04];
%    LonDMS=[-3 29  8.58;  7 48 4.44];
%    Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
%    Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
%    [x,y,utmzone] = deg2utm(Lat,Lon)
%
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06, Aug/06
% Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern 
%    hemisphere coordinates. 
% Aug/06: corrected m-Lint warnings
%-------------------------------------------------------------------------
% Argument checking
%
error(nargchk(2, 2, nargin));  %2 arguments required
n1=length(Lat);
n2=length(Lon);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end
% Memory pre-allocation
%
x=zeros(n1,1);
y=zeros(n1,1);
utmzone(n1,:)='60 X';
% Main Loop
%
for i=1:n1
   la=Lat(i);
   lo=Lon(i);
   sa = 6378137.000000 ; sb = 6356752.314245;
         
   %e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   %alpha = ( sa - sb ) / sa;             %f
   %ablandamiento = 1 / alpha;   % 1/f
   lat = la * ( pi / 180 );
   lon = lo * ( pi / 180 );
   Huso = fix( ( lo / 6 ) + 31);
   S = ( ( Huso * 6 ) - 183 );
   deltaS = lon - ( S * ( pi / 180 ) );
   if (la<-72), Letra='C';
   elseif (la<-64), Letra='D';
   elseif (la<-56), Letra='E';
   elseif (la<-48), Letra='F';
   elseif (la<-40), Letra='G';
   elseif (la<-32), Letra='H';
   elseif (la<-24), Letra='J';
   elseif (la<-16), Letra='K';
   elseif (la<-8), Letra='L';
   elseif (la<0), Letra='M';
   elseif (la<8), Letra='N';
   elseif (la<16), Letra='P';
   elseif (la<24), Letra='Q';
   elseif (la<32), Letra='R';
   elseif (la<40), Letra='S';
   elseif (la<48), Letra='T';
   elseif (la<56), Letra='U';
   elseif (la<64), Letra='V';
   elseif (la<72), Letra='W';
   else Letra='X';
   end
   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log( ( 1 +  a) / ( 1 - a ) );
   nu = atan( tan(lat) / cos(deltaS) ) - lat;
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   ta = ( e2cuadrada / 2 ) * epsilon ^ 2 * ( cos(lat) ) ^ 2;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000;
   yy = nu * v * ( 1 + ta ) + Bm;
   if (yy<0)
       yy=9999999+yy;
   end
   x(i)=xx;
   y(i)=yy;
   utmzone(i,:)=sprintf('%02d %c',Huso,Letra);
end
end




















