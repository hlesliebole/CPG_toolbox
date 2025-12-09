function [tidetime,tidehgt]=getztide(StartDatenum,EndDatenum,TimeZone,Datum)

% for all options see https://api.tidesandcurrents.noaa.gov/api/prod/

% scripps pier tide station number
st='9410230';

% Datum options are 'msl','navd','mllw' 
% TimeZone options are 'gmt','lst','lst_ldt' ...

api='https://tidesandcurrents.noaa.gov/api/datagetter?';

% get tide time series starting year hour zero
%bd=datestr([yr,1,1,0,0,0],'yyyymmdd HH:MM');
bd=datestr(StartDatenum,'yyyymmdd HH:MM');
% end 23:00 on 12/31
%ed=datestr([yr,12,31,23,0,0],'yyyymmdd HH:MM');
ed=datestr(EndDatenum,'yyyymmdd HH:MM');

url = [api 'begin_date=' bd '&end_date=' ed...
    '&station=' st '&product=Hourly_height&datum=' Datum '&units=metric'...
    '&time_zone=gmt&interval=h&application=web_services&format=csv'];
options = weboptions('Timeout', 30);

% Get tide time series prediction using webread.
% Requesting data in format=csv format results in a matlab table object
% with local 'DateTime' strings in column 1 and (feet, NAVD) 'Prediction' 
% tide levels in column 2.
T = webread(url,options);

% convert the tide table object to a matlab time series object
%tide=timeseries(T.Prediction,datenum(T.DateTime));

% add to single time hgt vectors
tidetime=T.DateTime;
%tidehgt=vertcat(tidehgt,T.Prediction);
tidehgt=T.WaterLevel;

end
