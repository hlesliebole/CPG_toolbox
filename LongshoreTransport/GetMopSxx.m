function [wavetime,Sxx]=GetMopSxx(stn)

urlbase=...  % nearshore station
'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';

% the entire mop hindcast data set is a combination of the mop hindcast.nc
%  and the mop nowcast.nc file, so download both and cat together

% read hindcast file dates and wave height info
urlend = '_hindcast.nc'; 
dsurl=strcat(urlbase,stn,urlend); % full url path and filename
wavetime=ncread(dsurl,'waveTime'); % read wave times
% convert mop unix time to UTC
hctime=datetime(wavetime,'ConvertFrom','posixTime');
hcSxx=ncread(dsurl,'waveSxx'); % read wave heights

% now read nowcast file
urlend = '_nowcast.nc';
dsurl=strcat(urlbase,stn,urlend);% full url path and filename

% read in netcdf wave forecast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to UTC
wavetime=datetime(wavetime,'ConvertFrom','posixTime');
wavetime=vertcat(hctime,wavetime); % combine with hindcast times

% read nowcast data
waveSxx=ncread(dsurl,'waveSxx'); 
Sxx=vertcat(hcSxx,waveSxx);

end