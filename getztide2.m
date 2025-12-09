function [tidetime,tidehgt]=getztide(StartDatenum,EndDatenum,TimeZone,Datum,Product)

% for all options see https://api.tidesandcurrents.noaa.gov/api/prod/

% scripps pier tide station number
st='9410230';

% Datum options are 'msl','navd','mllw' 
% TimeZone options are 'gmt','lst','lst_ldt' ...
% Product options 'Hourly_height'  or 'predictions'

api='https://tidesandcurrents.noaa.gov/api/datagetter?';

% get tide time series starting year hour zero
%bd=datestr([yr,1,1,0,0,0],'yyyymmdd HH:MM');
bd=datestr(StartDatenum,'yyyymmdd HH:MM');
% end 23:00 on 12/31
%ed=datestr([yr,12,31,23,0,0],'yyyymmdd HH:MM');
ed=datestr(EndDatenum,'yyyymmdd HH:MM');

url = [api 'begin_date=' bd '&end_date=' ed...
    '&station=' st '&product=' Product '&datum=' Datum '&units=metric'...
    '&time_zone=gmt&interval=h&application=web_services&format=csv'];
options = weboptions('Timeout', 60);

% Get tide time series prediction using webread.
% Requesting data in format=csv format results in a matlab table object
% with local 'DateTime' strings in column 1 and (feet, NAVD) 'Prediction' 
% tide levels in column 2.
T = webread(url,options);

% convert the tide table object to a matlab time series object
%tide=timeseries(T.Prediction,datenum(T.DateTime));

% add to single time hgt vectors
if ismember('DateTime',T.Properties.VariableNames)
tidetime=T.DateTime;
if strcmp(Product,'predictions')    
tidehgt=T.Prediction;
else
tidehgt=T.WaterLevel;
end
else
    fprintf('Error Using NOAA tide API with getztide2.m\n')
    tidetime=[];
    tidehgt=[];
end

end
