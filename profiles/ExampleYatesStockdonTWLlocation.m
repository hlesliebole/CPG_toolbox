

% uses functions
%    SA2NearestPointProfiles.m
%    intersections.m

MopNumber=582;

% get the most recent 1m spatial average data for this Mop
%   in the struct array SA
eval(['load(' sprintf('"M%5.5iSA.mat"',MopNumber) ',''SA'')']);

% most recent survey datenum
SurveyDatenum=SA(end).Datenum;

% get the the most recent profile from SA using the nearest points method
[X1D,Z1D]=SA2NearestPointProfiles(SA(end),25,4);

% get the xshore location of the MHW contour
xMSL=intersections([X1D(1) X1D(end)],[0.774 0.774],X1D,Z1D);
if numel(xMSL) == 0
    xMSL=X1D(find(~isnan(Z1D),1,'last'));
else
    XMSL=xMSL(1);
end

%% get recent waves
stn =['D' num2str(MopNumber,'%4.4i')]; % mop namefor waves
urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';  
% Set 'base' THREDDS URL and pieces to concatenate with user-defined station number (set above)
urlend = '_nowcast.nc';
dsurl = strcat(urlbase,stn,urlend);
% read in netcdf wave hs forecast convert to feet
wavehs=ncread(dsurl,'waveHs');
E=(wavehs/4).^2; % convert to energy
% read in netcdf wave forecast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime');
% reduce to hours since the last survey (assumed to be 0400 UTC on survey
% date)
surveytime=datetime(SurveyDatenum+4/24,'convertfrom','datenum');
idx=find( wavetime >= surveytime);
wavetime=wavetime(idx); % hourly datetime since last survey
E=E(idx); % hourly wave energy since last survey

%% run yates using MHW to shift entire profile in xshore direction
%  (aka hourly xshift time series)

%  Y09 Model

% Model Settings from Y09 paper Table 3 for Torrey Pines. Not sure
%  if these are precisely what is used in Fig 4.
b=0.07; % don't think "b" is explicitly given in paper. Estimated
        % from x=0 crossing of best fit line in Figure 3.
a=-0.0045; % from Table 3
Cminus=-1.38; % from Table 3
Cplus=-1.16; % from Table 3

% preallocate array large enough to hold results
Xshift=zeros(1,length(E)); % MSL shoreline position array

% set the first shoreline location as xMHW of survey
%  Y09 Fig 4 model starting on first data point
Xshift(1)=xMSL;

% step through hourly wave data and change the shoreline location

for i=2:numel(E)
    Eeq=a*Xshift(i-1)+b; % Yates eq. 4 equilibrium E
    deltaE=E(i)-Eeq; % Yates eq. 3 current diff from equilibrium E
    if( deltaE > 0)   
        dSdt=Cminus*sqrt(E(i))*deltaE; % Yates eq. 2 erosion
    else     
        dSdt=Cplus*sqrt(E(i))*deltaE; % Yates eq. 2 accretion
    end
    
    %Xshift(i)=Xshift(i-1)+dSdt; % change location
    Xshift(i)=dSdt;
end

Xshift=-Xshift;

%% make Stockdon tide + R2% TWL time series


%% step hourly through time since last survey, shift profile in
%   and find TWL x-intersection with the profile. If TWL exceeds 
%   max profile elevation, set water line location to x=0.



% get waves from the last survey to now
datestr(SA(end).Datenum)

