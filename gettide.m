%
%  CCCIA Imperial Beach "Dry Beach Forecast"
%
%  Beta code to make a realtime short-range forecast of the dry beach width
%  north of imperial beach pier.
%    - grabs latest tide info from NOAA
%    - grabs latest MOP Hs forecast from CDIP
%    - guesstimates tide+runup total water levels
%    - guesstimates resulting dry beach widths (if any)
%       based on a simplified two-slope profile assumption for IB beaches.
%    - makes user-friendly dry beach forecast product for stakeholders
%      (drybeach.png)
% 
%----------------------------------------------------------------------
% Step 1: Get the latest tide prediction from noaa
%----------------------------------------------------------------------

api='https://tidesandcurrents.noaa.gov/api/datagetter?';
% get tide time series starting -2 days from now
bd=datestr(now-2,'yyyymmdd HH:MM');
% and ending 3 days from now. Can be extended as desired.
ed=datestr(now+3,'yyyymmdd HH:MM');
% scripps pier tide station number
st='9410230';

% Build url string: note datum, units and time zone options.
% For all tide api options see https://tidesandcurrents.noaa.gov/api/
% Note: I'm just using the predcition here, so no anomaly info.
% This could be refined to grab both their actual water level product
% and the predictions to include the anomaly in the forecast.


url = [api 'begin_date=' bd '&end_date=' ed...
    '&station=' st '&product=predictions&datum=navd&units=english'...
    '&time_zone=lst_ldt&interval=h&application=web_services&format=csv'];
options = weboptions('Timeout', 30);

% Get tide time series prediction using webread.
% Requesting data in format=csv format results in a matlab table object
% with local 'DateTime' strings in column 1 and (feet, NAVD) 'Prediction' 
% tide levels in column 2.
T = webread(url,options);
% convert the tide table object to a matlab time series object
tide=timeseries(T.Prediction,datenum(T.DateTime));

%----------------------------------------------------------------------
%  Step 2:  get the latest MOP wave height forecast
%----------------------------------------------------------------------

% pick a mop number
stn='D0054';

% create url pathname to CDIP thredds netcdf file
urlbase=...
 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';
urlend='_forecast.nc';
% put the url string together
dsurl = strcat(urlbase,stn,urlend);
% read in netcdf wave hs forecast convert to feet
wavehs=ncread(dsurl,'waveHs')*3.28;
% read in netcdf wave forecast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime',...
    'TimeZone','America/Los_Angeles');

% combined mop wave hs and time vectors into a matlab timeseries object
wave=timeseries(wavehs,datenum(wavetime));

%----------------------------------------------------------------------
% Step 3: synchronize the tide and wave time series using the matlab
%         time series synchronize command with the "union" option.
%         In the future, dynamic profile variable time series (eg. msl 
%         shoreline migration based on Yates) can be similarly synced with 
%         the tide and wave times series for use in the final dry 
%         beach width calcs.
%----------------------------------------------------------------------

[tide,wave]=synchronize(tide,wave,'union');

%----------------------------------------------------------------------
% Step 4: Calculate a total water level based on tide and wave runup
%----------------------------------------------------------------------

%  For now, total water level "model" is just the sum of the tide elev 
%  and Hs at 10m depth.  This can be replaced by Stockdon runup
%  or an adjusted wave elevation based on historical flooding events etc.
%
%  If you want to make fake inundation event plots to see what it
%  might look like when we get storms this winter, you can multiply
%  wave.data timeseries below by a factor of 3 or 4.

totalH2O=tide.data+wave.data;

%----------------------------------------------------------------------
% Step 5: A simple I.B. "two-slope" dry beach width model to get started
%----------------------------------------------------------------------

%  Simple IB "two-slope" profile assumptions
%   1. I.B. profile can be approximated as two "fixed" (rarely need
%      to be modified) slopes, where the steeper beach face slope can
%      migrate landward/seaward over time. The two slopes are:
%       - a ~flat usually-dry upper beach extending from the houses/backbeach 
%         to the crown at the top of the usually-wet beach face.
%       - a steeper usually-wet beach face extending from the crown
%         seaward into the always-wet nearshore zone.
%   2.  While the usually-dry beach is ~flat or even slightly lower
%         than the crown, it is assigned a mild "effective slope" to
%         account for sand porosity when wave runup goes over the crown.
%         This adds a quasi-dynamic inundation aspect to the model vs
%         a bathtub prediction that would instantly reach the plaza
%         with even the smallest total water level elevation exceedance
%         of the crown elevation.
%
% Fixed model parameters:
%  1. mean slope of usually-wet beach face
%  2. "effective slope" of the usually-dry upper beach
%
% Dynamic profile parameters (updated by surveys or Yates predictions):
%  1. elevation of usually-dry, ~flat upper beach
%  2. xshore location of crown (eg. defined by subaerial survey) 
%     OR msl shoreline location on beach face (eg. defined by a survey or
%     estimated by Yates model).  These are *not* independent variables
%     owing to our use of a fixed beach face slope, so you can use
%     one or the other, or both with the beach slope value to work
%     out best fit values to use when making a two-slope profile update.

% --- start initial 4 profile settings ----

% We need twoslope settings. Made these up as a starting point.  A better 
% mean beach face slope can be derived from surveys.  
% The effective back beach slope should be the last thing to be tuned 
% based on historical back beach and plaza wave inundation events.
BeachFaceSlope=1/30;
EffectiveBackSlope=1/100;

%  Next, we need a current mean elevation (feet, NAVD) of the 
%  usually-dry beach (same as height of crown in our simple 2-slope profile 
%  shape). Guesstimated typical upper beach height from moff and nichols 
%  nourishment planning report, recent surveys can refine.

BeachHeight=12;

% Finally we need a distance in feet from houses/seawall/beginning of sand 
% to either the crown of beach (the slope break in our two-slope profile)
% OR the msl shoreline location of the beach face.  I'm going with a
% beach crown location for starters (guesstimated this from current 
% surfline IB web cam images and the handy distance scale on google maps 
% as a ballpark starting point. Should be updated going forward based on 
% latest survey info)

CrownDistance=200;

% --- end of profile settings ----

% Now, based on crown location, get distance from the houses to the 
% msl shoreline location on the beach face. (Based on msl=2.92ft NAVD)

MSLshoreDistance=CrownDistance+(BeachHeight-2.92)/BeachFaceSlope;

% Note: for the profile to "evolve" between surveys based on mop
%  waves, we could use Yates to predict a shift in this msl shoreline
%  location, which in turn can be used to calculate a new crown
%  location (CrownDistance) based on the fixed BeachFaceSlope.

%----------------------------------------------------------------------
% Step 6: Calculate a dry beach width time sereis based on two-slope 
%         profile and total water level predictions
%----------------------------------------------------------------------

% if total water level is less than BeachHeight, then the swash is
% restricted to the beach face.
i=find(totalH2O <= BeachHeight);
DryWidth(i)=MSLshoreDistance-((totalH2O(i)-2.92)/BeachFaceSlope);
% if total water level is above than BeachHeight then water
%  extends up the effective backbeach slope. Negtive dry beach
%  values (flooding of area behind beach, aka onto the plaza) is possible 
%  with large enough total water level.
j=find(totalH2O > BeachHeight);
DryWidth(j)=CrownDistance-((totalH2O(j)-BeachHeight)/EffectiveBackSlope);

%----------------------------------------------------------------------
%  Step 7:  Make A User-Friendly Plot of the Dry Beach Width forecast
%----------------------------------------------------------------------

figure('position',[450   350   850   425]);
ymin=-100;ymax=800;
xmin=tide.time(1);xmax=tide.time(end);
% highlight beach and beach face xshore locations 
fill([xmin xmax xmax xmin],[0 0  CrownDistance CrownDistance],...
    'y','edgecolor','none');
hold on;
fill([xmin xmax xmax xmin],[CrownDistance CrownDistance ymax ymax],...
    [.7 .7 0],'edgecolor','none');
% highlight plaza xshore location
fill([xmin xmax xmax xmin],[0 0 ymin ymin],...
    'r','edgecolor','none');
% plot total water level converted to dry beach width
plot(tide.time,DryWidth,'b-','linewidth',2);
fill([xmin tide.time' xmax],[ymax DryWidth ymax],'b','edgecolor','none');
% show msl shoreline position
plot([xmin xmax],[MSLshoreDistance MSLshoreDistance],'k--');

% label xaxis with days
datetick('x','ddd mm/dd');
% flip yaxis upside down and set fix axes ranges
set(gca,'ydir','reverse','xlim',[xmin xmax],'ylim',[ymin ymax]);

% label things 
title('Imperial Beach: Experimental Dry Beach Width Forecast',...
    'fontsize',14,'fontweight','bold');
xlabel('Local Date/Time','fontsize',12,'fontweight','bold');
ylabel('Dry Beach Width (ft)','fontsize',12,'fontweight','bold');
text(xmin+.1, CrownDistance,'Beach Crown','verticalalign','top');
text(xmin+.1, MSLshoreDistance,'MSL Shoreline','verticalalign','top');
text(xmin+.6, CrownDistance/2,'Usually Dry Beach','verticalalign','bottom');
text(xmin+.6, CrownDistance+(MSLshoreDistance-CrownDistance)/2,...
    'Beach Face','verticalalign','bottom');
text(xmin+.6, ymax-100,'Tide + Wave Runup Wet Beach','verticalalign','bottom');
hold on;
text(xmin+.1, 0,'Palm Plaza','verticalalign','bottom'); 
% mark the current time with a vertical line
plot([datenum(now) datenum(now)],[ymin ymax],'k-');
text(datenum(now),ymax-20,datestr(now,' HH AM'),'verticalalign','bottom',...
    'horizontalalign','left','fontsize',12,'fontweight','bold')

%----------------------------------------------------------------------
% Finally, make a png image that could be used on the website
%----------------------------------------------------------------------

pfile=['beachwidth.png'];
print(gcf,'-dpng','-r300','-loose',pfile);