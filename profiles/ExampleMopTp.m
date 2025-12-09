%% example code to get latest Tp at Fletcher

stn = 'D0654'; % fletcher

% contruct name of mop realtime file
urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';  
% Set 'base' THREDDS URL and pieces to concatenate with user-defined station number (set above)
urlend = '_nowcast.nc'; % change this to _forecast.nc to make forecast prediction
ncfile=strcat(urlbase,stn,urlend);

% option to display all the variables in the file to teh screen
%ncdisp(ncfile)

% get the mop data times
wavetime=ncread(ncfile,'waveTime');
% convert time to a matlab datetime
wavetime=datetime(wavetime,'ConvertFrom','posixtime');

% get Mop Tp values
Tp=ncread(ncfile,'waveTp');

% print the most recent time and value
fprintf('\nLatest Mop Tp: %s UTC   %f4.1 s\n',wavetime(end),Tp(end))

