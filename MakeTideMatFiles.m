%
%  Makes LaJollaMslTide20002022.mat
%
%   containing the stuct array Tide with the fields
%      Tide.Datenum, Tide.MSL
% 

api='https://tidesandcurrents.noaa.gov/api/datagetter?';

% scripps pier tide station number
st='9410230';

% Build url string: note datum, units and time zone options.
% For all tide api options see https://tidesandcurrents.noaa.gov/api/
% Note: I'm just using the predcition here, so no anomaly info.
% This could be refined to grab both their actual water level product
% and the predictions to include the anomaly in the forecast.

Tide.Datenum=[];
Tide.MSL=[];

for yr=2000:2022
    
fprintf('%i\n',yr)

% get tide time series starting year hour zero
bd=datestr([yr,1,1,0,0,0],'yyyymmdd HH:MM');
% end 23:00 on 12/31
ed=datestr([yr,12,31,23,0,0],'yyyymmdd HH:MM');

url = [api 'begin_date=' bd '&end_date=' ed...
    '&station=' st '&product=predictions&datum=msl&units=metric'...
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
Tide.Datenum=vertcat([Tide.Datenum],T.DateTime);
Tide.MSL=vertcat([Tide.MSL],T.Prediction);

end

% figure;
% plot(Tide.Datenum,Tide.MSL,'b-','linewidth',1);
% datetick('x','mm/yy');

save LaJollaMslTide20002022.mat Tide

