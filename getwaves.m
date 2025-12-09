function [wavetime,wavehs]=getwaves(MopNumber,StartDatenum,EndDatenum,TimeZone)

if strcmpi(TimeZone,'lst') 
    tz='America/Los_Angeles';
else
    tz='Etc/UTC';
end

%fprintf('Getting Wave Data for Mop %i of 664...\n',mop)
stn=['D0' num2str(MopNumber)];

urlbase=...  % nearshore station
'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';

% the entire mop hindcast data set is a combination of the mop hindcast.nc
%  and the mop nowcast.nc file, so download both and cat together

% read hindcast file dates and wave height info
urlend = '_hindcast.nc'; 
dsurl=strcat(urlbase,stn,urlend); % full url path and filename
wavetime=ncread(dsurl,'waveTime'); % read wave times
% convert mop unix time to UTC
hctime=datetime(wavetime,'ConvertFrom','posixTime','TimeZone',tz);
hchs=ncread(dsurl,'waveHs'); % read wave heights
% hctp=ncread(dsurl,'waveTp');
% hcta=ncread(dsurl,'waveTa');
% hcdp=ncread(dsurl,'waveDp');
% hcSxx=ncread(dsurl,'waveSxx');
% hcSxy=ncread(dsurl,'waveSxy');

% now read nowcast file
urlend = '_nowcast.nc';
dsurl=strcat(urlbase,stn,urlend);% full url path and filename

% read in netcdf wave forecast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to UTC
wavetime=datetime(wavetime,'ConvertFrom','posixTime','TimeZone',tz);
wavetime=vertcat(hctime,wavetime); % combine with hindcast times

% read nowcast data
wavehs=ncread(dsurl,'waveHs'); 
% wavetp=ncread(dsurl,'waveTp');
% waveta=ncread(dsurl,'waveTa');
% wavedp=ncread(dsurl,'waveDp');
% waveSxx=ncread(dsurl,'waveSxx');
% waveSxy=ncread(dsurl,'waveSxy');


% combine with hindcast file
wavehs=vertcat(hchs,wavehs); 
% wavetp=vertcat(hctp,wavetp);
% waveta=vertcat(hcta,waveta);
% wavedp=vertcat(hcdp,wavedp);
% waveSxx=vertcat(hcSxx,waveSxx);
% waveSxy=vertcat(hcSxy,waveSxy);

%bchnorm=ncread(dsurl,'metaShoreNormal');

idx=find(datenum(wavetime) >= StartDatenum & datenum(wavetime) <= EndDatenum);

wavetime=wavetime(idx);
wavehs=wavehs(idx);

end

