function [wavetime,Sxx,Sxy]=GetCpgMopRadiationStresses(MopID,FreqStart,FreqEnd,NormalRotation)

% Reads Mop spectral data from CDIP online netcdf files and calculates
%  Sxy over a specified freq range with an optional rotation of the official
%  shore normal

% Input:
%
% MopID = Can be either CDIP county-number character string (eg. 'D0654') 
%           or its equivalent sequential CA numeric number (eg. 654) 
% 
% FreqStart,FreqEnd = Starting and ending frequency bands to include 
%                  (Mop spectra have a frequency range of 0.04 - 0.40Hz)
% 
% NormalRotation= desired rotation of the official shore normal in degrees 
%         (eg. -5 or 5) in the true compass heading coordinate frame, so a 
%          positive rotation will turn a west facing beach normal more northward.
% 
% Mop spectra have a frequency range of 0.04 - 0.40Hz, so
% 
% [wavetime,Sxx,Sxy]=GetCpgMopSxy(MopID,0.04,0.40,0.0) 
% 
% should return Sxy values matching the netcdf "waveSxx" and "waveSxy" parameters.

% Output:
%
% wavetime = Matlab datetimes of the hourly Sxy data
%
% Sxy = Sxy values (m^2)

%% Convert MopID to a CDIP county-number character string
%   to access netcdf file. If MopID is already a county-number
%   character string it is just returned as MopName

[MopName,MopNumber]=GetMopNameNumber(MopID);

%% ---------------------------------------------------------------
%% Get spectral data from CDIP nedcdf hindcast and nowcast files
%% ---------------------------------------------------------------

urlbase=...  % nearshore station
'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';

% the entire mop hindcast data set is a combination of the mop hindcast.nc
%  and the mop nowcast.nc file, so download both and cat together

% read hindcast file dates and wave height info
urlend = '_hindcast.nc'; 
dsurl=strcat(urlbase,MopName,urlend); % full url path and filename
hctime=ncread(dsurl,'waveTime'); % read wave times
hctime=datetime(hctime,'ConvertFrom','posixTime');% convert mop unix time to UTC
hca0=ncread(dsurl,'waveEnergyDensity'); % a0 wave energy
hca2=ncread(dsurl,'waveA2Value'); % a2 dir fourier coeff
hcb2=ncread(dsurl,'waveB2Value'); % b2 dir fourier coeff

% now read nowcast file data
urlend = '_nowcast.nc';
dsurl=strcat(urlbase,MopName,urlend);
wavetime=ncread(dsurl,'waveTime');
wavetime=datetime(wavetime,'ConvertFrom','posixTime');
wavea0=ncread(dsurl,'waveEnergyDensity');
wavea2=ncread(dsurl,'waveA2Value');
waveb2=ncread(dsurl,'waveB2Value'); 

% combine the hindcast and nowcast data
wavetime=vertcat(hctime,wavetime); % combine with hindcast times
e=vertcat(hca0',wavea0')';
a2=vertcat(hca2',wavea2')';
b2=vertcat(hcb2',waveb2')';

%% get the freq band info and offical Mop shore normal
%   from the nowcast netcdf file

freq=ncread(dsurl,'waveFrequency');
bw=ncread(dsurl,'waveBandwidth');
snorm=ncread(dsurl,'metaShoreNormal');

%% rotate the normal by the input value
snorm=snorm+NormalRotation;

%% ---------------------------------------------------------------
%% Calculate Sxx and Sxy
%% ---------------------------------------------------------------

% calculate C and Cg for each frequency band 
c=NaN(1,20);Cg(1,20)=NaN; % preallocate

d=10.0; % mop water depth
for n=1:length(freq) % loop through freq bands
 f=freq(n);
 % C.S. Wu linear dispersion algorithm to get radian wavenumber
 %  based on freq and depth
 a=4.0243*f*f*d;
 yhat=a.*(1+1.26*exp((-1.84).*a));t=exp((-2).*yhat);
 aka=a;aka(a >= 1)=a(a >= 1).*(1+2*t(a >= 1).*(1+t(a >= 1)));
 aka(a < 1)=sqrt(a(a < 1)).*(1+(a(a < 1)./6).*(1+a(a < 1)./5));
 k=abs(aka./d);
 % phase and group velocities
 c(n)=f*2*pi/k;
 Cg(n)=pi*(f./k).*(1+2*k.*d./sinh(2*k.*d));
end

% calculate n=Cg/C vs freq for Sxx and Sxy equations  
nC=Cg./c;  

% rotate dir fourier coeffs relative to the shore normal 
%  using trig rotation formulas
s=sind(-snorm);
c=cosd(-snorm);
s2=sind(-2*snorm);
c2=cosd(-2*snorm);
%b1r=b1*c+a1*s;
%a1r=a1*c-b1*s;
b2r=b2*c2+a2*s2;
a2r=a2*c2-b2*s2;

% figure out the frequency range to include based on input starting and 
% ending angles.  Find the Mop freq band index closest to each
[~,ifmin]=min(abs(freq-FreqStart)); % FreqStart Mop freq index
[~,ifmax]=min(abs(freq-FreqEnd));% FreqEnd Mop freq index

nf=ifmin:ifmax; % Mop freq range indexes

% Calculate Sxy from rotated b2, need to flip sign for minus=southward
%  transport convention
Sxy=-0.5*sum(nC(nf)'.*e(nf,:).*bw(nf).*b2r(nf,:),'omitnan')';

% Calculate Sxx from rotated a2
nCa=repmat(nC',1,size(a2r,2));
Sxx=sum(e(nf,:).*bw(nf).*(nCa(nf,:)-0.5+nC(nf)'.*(0.5*(1+a2r(nf,:)))),'omitnan')';

end

%% supporting functions

function [MopName,MopNumber]=GetMopNameNumber(MopID)

% A function to make it easy to work with both Mop character names
% based on the County, or their state-wide sequential Mop number. 
% Whichever is input to the function, it returns it and its name/number 
% counterpart as MopName and MopNumber.

County=['D    ';'OC   ';'L    ';'VE   ';'B    ';'SL   ';'MO   ';'SC   ';...
             'SM   ';'SF   ';'MA   ';'SN   ';'M    ';'HU   ';'DN   '];
CountyBeginNumber=[0 1210 1878 3063 3741 5529 6281 7202 ...
                  7530 7977 8066 8659 9172 10225 11221];

if isnumeric(MopID) % figure out MopName from the input MopNumber

    MopNumber=MopID; 
    % figure out what county the state number is in and its county number
    CountyNumbers=MopID-CountyBeginNumber;
    CountyNumbers(CountyNumbers < 1)=NaN;
    [CountyNumber,CountyIndex]=min(CountyNumbers); 
    % build MopName from its county character(s) and number
    MopName=County(CountyIndex,:);
    CountyNumeric=numel(find(MopName == ' '));
    if numel(find(MopName == ' ')) == 3
      MopName(3:5)=num2str(CountyNumber,'%3.3i');
    else
      MopName(2:5)=num2str(CountyNumber,'%4.4i');
    end

else % figure out MopNumber from the input MopName

    MopName=MopID;
    % find county that matches characters in MopName
    CountyIndex=find(strcmp(cellstr(County(:,1:5)),...
        cellstr(regexprep(MopName,'[0-9]',' '))));
    % extract the county number from MopName
    CountyNumber=str2num(regexprep(MopName,'[A-Z]',' '));
    % convert county number to state number
    MopNumber=CountyBeginNumber(CountyIndex)+CountyNumber;
    
end

end
