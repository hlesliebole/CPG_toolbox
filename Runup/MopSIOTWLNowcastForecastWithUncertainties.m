%% Makes Total Water Level Tide + standard Stockdon R2% runup estimates 
% from a MOP frequency spectrum time series and NOAA tide stations.
%
%  This version loads tide predictions from 2000 through 2025 in the
%  file CAtidePrediction.mat so it can be used for a few years.
%

% -Uses seasonal Mop beach slope sine wave from Susheel's analysis saved 
%   in the table "VosT" in the file VosTransects-filt.mat

% -Gets Mop predictions from the CDIP thredds server. Unshoals the
%  frequency spectra to deep water and calculates Ho, Tp and Lo @ Tp.

% -Plots runup uncertainty with the slope, using the min-max
%   possible slopes of the seasonal slope sine wave. 

% -Plots runup uncertainty with Tp (aka deep water wavelength Lo), using
%  the peak period of the spectrum in the swell band (T > 10s), and the
%  sea band (T < 10s ). If either band does not have an actual peak its
%  Tp value will be at the T = 10s sea-swell boundary.

%% Settings
% MopID can be County-number character string in quotes or a state-wide sequential 
% mop number not in quotes.
%MopID='OC399'; % Balboa Pier
% MopID='D0662'; % N Cardiff parking lot
MopID='D0511'; % Fletcher Cove, 200m north of northward looking cam
% MopID='D0640'; % dog beach

% Rstart=datenum(2024,8,1);  % runup time series start date
% Rend=datenum(2024,9,10); % runup time series end date
%% --------------------------------------------------------------
%% Set time range of prediction: most recent survey date (@4am UTC) to now
Rstart=floor(datenum(datetime('now','timezone','America/Los_Angeles')))-14; 
% end a day into the future. Could go out 7 days or so with full forecast
Rend=floor(datenum(datetime('now','timezone','America/Los_Angeles')))+7; 
% Rstart=datenum(2022,10,1);  % runup time series start date
% Rend=datenum(2023,4,1); % runup time series end date
MopQCflag='on'; 
%  MopQCflag ['on'/'off'] quality control flag; 
%  When set 'on' Ho = 0 will be returned for times when either the 
%  available buoys for the swell band or the sea band predictions were 
%  deemed to be too far away compared to the planned buoy data coverage 
%  for this Mop location.  Set to 'off' it will return results no matter
%  what most of the time.


%% 1. load Susheel's table variable with seasonal slope coefficients 
load VosTransects-filt.mat

%% 2. Load CA tide info. Tide heights are MSL, tide times are UTC
%   has 4 variables 1) CAtidePrediction(8x227928) hourly tides predictions for
%   2000-2025 (predictions, no anomalies) from 8 tide stations along
%   coast. 2) MopTideStations(1x8) Mop numbers associated with tide
%   stations. 3) TideDatetimeGMT(227928x1) hourly tide datetimes, and
%   4) Tidestations(1x8) NOAA tide station numbers.  SF, Humboldt and
%   Crescent City harbor tide info was not used because the alongcoast 
%   phase lags and amplitudes seemed a little off compared to the 
%   other 8 more exposed stations.
load CAtidePrediction.mat

% convert the MopID to a CA sequential Mop number to find slope parameters
[MopName,MopNumber]=GetMopNameNumber(MopID);

%% 3. Slope Parameters
%   Either define a fixed slope or try and get seasonal slope
%   function parameters based on Susheel's satellite analysis.

A=NaN;B=NaN;S=NaN; % default to get Susheel slopes 

% If you want to overide the Susheel slopes with a fixed slope, set B to
%    the slope you want and A and S to zero. Otherwise comment out 
%    this fixed slope line 
A=0;B=0.15;S=0; % example for fixed Solana Pad Beach Face slope =0.15

%% ------------------------

% Use Susheel slopes if B is a NaN
if isnan(B)

% find Susheel slope params that are closest to this Mop number
[mopdist,idx]=min(abs([VosT.BackMop{:}]- MopNumber));
else
    mopdist=0;
end

% if params are within 3 Mop number "lengths" of the Mop, proceed
if mopdist <= 3 

  if isnan(B)
% slope seasonal sine wave parameters for Beta= B + A*cos(2*pi*(doy'+S)/365); 
   idx=find(round([VosT.BackMop{:}]) == MopNumber);
   A=VosT.A(idx);
   B=VosT.B(idx); 
   S=VosT.S(idx);
  end

  % if there is valid slope information at this Mop, proceed
  if ~isnan(S) 

%% Get MOP Ho,Lo predictions for runup date range
[moptime,fidx,Ho,Lo,LoSwl,LoSea,Tp,TpSwl,TpSea]=GetMopStockdonParams(MopName,Rstart,Rend,MopQCflag);

%% Use slope parameters with the day of the year to calculate
%   the seasonal slope time series to go with the Ho Lo time series
doy=day(datetime(moptime,'convertfrom','datenum'),'dayofyear');
Beta= B + A*cos(2*pi*(doy'+S)/365); % Seasonal slope sine wave time series
BetaMin=(B-A)*ones(size(moptime))'; % Min satellite slope at all times
BetaMax=(B+A)*ones(size(moptime))'; % Max satellite slope at all times

%% Stockdon Eq 19 runup elevation time series
R2=1.1*(0.35*Beta.*sqrt(Ho.*Lo)+0.5*sqrt(Ho.*Lo.*(0.563*(Beta.^2)+0.004)));

% R2 min-max with min-max seasonal slopes
R2min=1.1*(0.35*BetaMin.*sqrt(Ho.*Lo)+0.5*sqrt(Ho.*Lo.*(0.563*(BetaMin.^2)+0.004)));
R2max=1.1*(0.35*BetaMax.*sqrt(Ho.*Lo)+0.5*sqrt(Ho.*Lo.*(0.563*(BetaMax.^2)+0.004)));
% R2 min-max using sea-swell Tp [Lo] values 
R2swl=1.1*(0.35*Beta.*sqrt(Ho.*LoSwl)+0.5*sqrt(Ho.*LoSwl.*(0.563*(Beta.^2)+0.004)));
R2sea=1.1*(0.35*Beta.*sqrt(Ho.*LoSea)+0.5*sqrt(Ho.*LoSea.*(0.563*(Beta.^2)+0.004)));

%% Get matching Tide time series for this MopNumber and date range
%   tide elevations are linearly interpolated between tide stations based
%   on the mop number.  If south of SIO pier or north of Arena Cove,
%   linear extrapolatation is used.
if(MopNumber >= min(MopTideStations) && MopNumber <= max(MopTideStations))
  MopTide=interp2(datenum(TideDatetimeGMT),MopTideStations,CAtidePrediction,moptime,MopNumber);
else
   % need to extrapolate
  idx=find(ismember(datenum(TideDatetimeGMT),moptime));
  for i=1:numel(idx)
    MopTide(i)=interp1(MopTideStations,CAtidePrediction(:,idx(i)),MopNumber,'linear','extrap');
  end
end

%% TWL time series
TWL=R2+MopTide;
TWLmin=R2min+MopTide;
TWLmax=R2max+MopTide;
TWLsea=R2sea+MopTide;
TWLswl=R2swl+MopTide;

%% Make some runup time series plots with uncertainties
%figure('position',[ 87          60        1186         737 ]);
figure('position',[ 87         252        1186         545 ]);

% % plot with slope uncertainty
% ax1=axes('position',[.05 .53 .9 .4]);
% p=plot([moptime'; moptime'],[TWLmin; TWLmax],'r-');
% hold on;
% idx=find(moptime <= now);
% r=plot(moptime(idx)',TWL(idx),'b.-');
% idx=find(moptime > now); % plot forecast part as red
% r2=plot(moptime(idx)',TWL(idx),'r.-');
% 
% datetick('x','mm/dd/yy','keepticks','keeplimits');grid on;
% set(gca,'xlim',[Rstart Rend]);
% set(gca,'fontsize',14);xlabel('Time (UTC)');ylabel('Tide + R2% (m, MSL)');
% if A == 0 && S == 0
%     txt=['Using Fixed Slope: ' num2str(B,'%5.3f')];
%     title({['Mop ' MopName ' Tide+Runup TWL Prediction with Fixed Slope'],txt},'fontsize',16);
%     legend(r,'TWL using Stockdon Runup with Seasonal Satellite Slope Estimate',...
%     'location','southoutside');
% else
%     txt=['Seasonal Slopes are: ' num2str(B,'%5.3f') ' ' char(177) ' ' num2str(A,'%5.3f')...
%     ' with max slope on ' datestr(datenum(2023,1,1)+S,'mmm-dd') ];
%     title({['Mop ' MopName ' Tide+Runup TWL Prediction with Slope Uncertainty'],txt},'fontsize',16);
%     legend([r p(1)],'TWL using Stockdon Runup with Seasonal Satellite Slope Estimate',...
%     'TWL uncertainty with Min-Max Range of Satellite Slopes at this Mop',...
%     'location','southoutside');
% end
% 
% 
% % plot with Lo uncertainty
% ax2=axes('position',[.05 .05 .9 .4]);

% switch to PST
moptime=moptime-8/24;
p=plot([moptime'; moptime'],[TWLsea; TWLswl],'-','color',[.7 .7 .7],'linewidth',2);
hold on;
%idx=find(moptime <= now);
r=plot(moptime',TWL,'b.-','linewidth',2);
%idx=find(moptime > now); % plot forecast part as red
r2=plot(moptime(fidx)',TWL(fidx),'r.-','linewidth',2);
tl=plot(moptime,moptime*0+2.1,'k--','linewidth',2);

set(gca ,'xtick',Rstart:Rend)
datetick('x','ddd mm/dd','keepticks','keeplimits');grid on;
set(gca,'xlim',[Rstart Rend],'ylim',[-1 5.5]);
set(gca,'fontsize',16);xlabel('Time (PST)');ylabel('Tide + R2% (m,MSL)');
% txt=['Seasonal Slopes are: ' num2str(B,'%5.3f') ' ' char(177) ' ' num2str(A,'%5.3f')...
%     ' with max slope on ' datestr(datenum(2023,1,1)+S,'mmm-dd') ];
txt=['Using Fixed Slope: ' num2str(B,'%5.3f')];
title({['Mop ' MopName ' Tide+ Stockdon Runup TWL Prediction with Tp [Lo] Uncertainty'],txt},'fontsize',16);
legend([r r2 p(1) tl],'TWL using Stockdon Runup with Fixed Slope 0.15',...
    'ECMWF-driven Forecast',...
    'TWL uncertainty using the Swell Tp [Lo] or Sea Tp [Lo]','Terrace Elevation',...
    'location','southoutside');
set(gca,'linewidth',2)
  
  else
    fprintf('No Satellite-derived beach slopes at this Mop location. %s\n',MopName)
  end
else
    fprintf('No Satellite-derived beach slopes found within 300m of Mop %s\n',MopName)
end

makepng('../../www/MopD0511RunupForecast.png')

% end of main script


% %% ---------------------------------------------------------
% 
% function [MopName,MopNumber]=GetMopNameNumber(MopID)
% 
% % A function to make it easy to work with both Mop character names
% % based on the County, or their state-wide sequential Mop number. 
% % Whichever is input to the function, it returns it and its name/number 
% % counterpart as MopName and MopNumber.
% 
% County=['D    ';'OC   ';'L    ';'VE   ';'B    ';'SL   ';'MO   ';'SC   ';...
%              'SM   ';'SF   ';'MA   ';'SN   ';'M    ';'HU   ';'DN   '];
% CountyBeginNumber=[0 1210 1878 3063 3741 5529 6281 7202 ...
%                   7530 7977 8066 8659 9172 10225 11221];
% 
% if isnumeric(MopID) % figure out MopName from the input MopNumber
% 
%     MopNumber=MopID; 
%     % figure out what county the state number is in and its county number
%     CountyNumbers=MopID-CountyBeginNumber;
%     CountyNumbers(CountyNumbers < 1)=NaN;
%     [CountyNumber,CountyIndex]=min(CountyNumbers); 
%     % build MopName from its county character(s) and number
%     MopName=County(CountyIndex,:);
%     CountyNumeric=numel(find(MopName == ' '));
%     if numel(find(MopName == ' ')) == 3
%       MopName(3:5)=num2str(CountyNumber,'%3.3i');
%     else
%       MopName(2:5)=num2str(CountyNumber,'%4.4i');
%     end
% 
% else % figure out MopNumber from the input MopName
% 
%     MopName=MopID;
%     % find county that matches characters in MopName
%     CountyIndex=find(strcmp(cellstr(County(:,1:5)),...
%         cellstr(regexprep(MopName,'[0-9]',' '))));
%     % extract the county number from MopName
%     CountyNumber=str2num(regexprep(MopName,'[A-Z]',' '));
%     % convert county number to state number
%     MopNumber=CountyBeginNumber(CountyIndex)+CountyNumber;
% 
% end
% 
% end
% 
% %% ---------------------------------------------------------
% 
% function [moptime,fidx,Ho,Lo,LoSwl,LoSea,Tp,TpSwl,TpSea]=GetMopStockdonParams(MopID,Rstart,Rend,MopQCflag)
% 
% % first get the Mop frequency spectra
% [freq,bw,depth,moptime,fidx,a0] = GetMOPfreqSpectra(MopID,Rstart,Rend,MopQCflag);
% 
% % get unshoaling coefficients and deep water wavelength as a function
% %  of the Mop frequencies and depth
% 
% % need group speed Cg at Mop depth
% [Lf,Cf,Cgf]=LinearDispersion(freq,depth);
% % need group speed Cgo and wavelengths Lo in deep water
% [Lof,Cof,Cgof]=LinearDispersion(freq,2000);
% 
% % unshoal Mop freq spectra to deep water and convert to m^2 energies
% a0deep=a0.*repmat((Cgf./Cgof),1,size(a0,2));
% 
% % get deep water significant wave height
% Ho=4*sqrt(sum(a0deep.*repmat(bw,1,size(a0,2)),'omitnan'));
% 
% % find peak energy density frequency band for each spectrum
% [Epeak,fpeak]=max(a0deep);
% [EpeakSwl,fpeakSwl]=max(a0deep(1:10,:));
% [EpeakSea,fpeakSea]=max(a0deep(11:end,:));fpeakSea=fpeakSea+10;
% 
% % use peak freq bands to define Tp and Lo for each spectrum
% T=1./freq;
% Tp=T(fpeak)';
% TpSwl=T(fpeakSwl)';
% TpSea=T(fpeakSea)';
% Lo=Lof(fpeak)';
% LoSwl=Lof(fpeakSwl)';
% LoSea=Lof(fpeakSea)';
% 
% end
% 
% %% ------------------------------------------------------------------
% function [freq,bw,depth,moptime,fidx,a0] = GetMOPfreqSpectra(MopName,dstart,dend,MopQCflag)
% 
% % returns spectral data from the nowcast and ecmwf forecast thredds
% %  netcdf files for the mop specified by MopID and dstart to dend datenumber time range
% 
% % create urls for MopID on thredds server
% urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';  
% dset=strcat(urlbase,MopName,'_nowcast.nc'); % nowcast url
% dset2=strcat(urlbase,MopName,'_ecmwf_fc.nc'); % forecast url
% % Display variable names in the netcdf file in the command window
% %ncdisp(dset,'/','min')
% 
% % load nowcast time variable, find time indices needed
%   wtime = ncread(dset,'waveTime');
%   % convert posix time to datenumber
%   wtime = double(wtime)/86400 + datenum(1970,1,1);
%   datestr(wtime(end))
% 
%   % load forecast time variable, find time indices needed
%   wtime2 = ncread(dset2,'waveTime');
%   % convert posix time to datenumber
%   wtime2 = double(wtime2)/86400 + datenum(1970,1,1);
% 
% %  option for datetime variable
% % wavedatetime=datetime(wtime,'ConvertFrom','posixTime');
% 
% %% get spectral data from the nowcast file
% 
% a0=[];a1=[];b1=[];a2=[];b2=[];moptime=[];
% % check if at least some of the data is in the nowcast file
% if dstart <= wtime(end) 
% 
%   dindices = find(wtime >= dstart & wtime <= dend);
%   dlength = size(dindices,1);
%   nstart_1d = dindices(1);
%   ncount_1d = dlength;
%   moptime = wtime(dindices(1):dindices(ncount_1d));
%   % load qc flags; records good where qcflag = 1 
%   qcflag = ncread(dset,'waveFlagPrimary',nstart_1d,ncount_1d);
%   mask_1d = qcflag==1;
% 
%    % load 2d var (a1), mask off before plotting
%   freq = ncread(dset,'waveFrequency');
%   bw = ncread(dset,'waveBandwidth');
%   depth=ncread(dset,'metaWaterDepth');
%   nstart_2d = [1, dindices(1)];
%   ncount_2d = [size(freq,1), dlength];
%   a0 = ncread(dset,'waveEnergyDensity',nstart_2d,ncount_2d);
%   % a1 = ncread(dset,'waveA1Value',nstart_2d,ncount_2d);
%   % a2 = ncread(dset,'waveA2Value',nstart_2d,ncount_2d);
%   % b1 = ncread(dset,'waveB1Value',nstart_2d,ncount_2d);
%   % b2 = ncread(dset,'waveB2Value',nstart_2d,ncount_2d);
%   mask_2d = repmat(mask_1d,1,size(freq,1))';
%   if strcmp(MopQCflag(1:2),'on')
%    a0(mask_2d~=1) = NaN;
%   end
%   % a1(mask_2d~=1) = NaN;
%   % b1(mask_2d~=1) = NaN;
%   % a2(mask_2d~=1) = NaN;
%   % b2(mask_2d~=1) = NaN;
% 
% end
% 
% %% get forecast file data
% 
% na0=[];na1=[];nb1=[];na2=[];nb2=[];nptime=[];
% 
% if dend >= wtime2(1)
% 
%   dindices = find(wtime2 >= dstart & wtime2 <= dend);
%   dlength = size(dindices,1);
%   nstart_1d = dindices(1);
%   ncount_1d = dlength;
%   nptime = wtime2(dindices(1):dindices(ncount_1d));
% 
%   % load qc flags; records good where qcflag = 1 
%   qcflag = ncread(dset,'waveFlagPrimary',nstart_1d,ncount_1d);
%   mask_1d = qcflag==1;
% 
%    % load 2d var (a1), mask off before plotting
%   freq = ncread(dset,'waveFrequency');
%   bw = ncread(dset,'waveBandwidth');
%   depth=ncread(dset,'metaWaterDepth');
%   nstart_2d = [1, dindices(1)];
%   ncount_2d = [size(freq,1), dlength];
%   na0 = ncread(dset,'waveEnergyDensity',nstart_2d,ncount_2d);
%   % na1 = ncread(dset,'waveA1Value',nstart_2d,ncount_2d);
%   % na2 = ncread(dset,'waveA2Value',nstart_2d,ncount_2d);
%   % nb1 = ncread(dset,'waveB1Value',nstart_2d,ncount_2d);
%   % nb2 = ncread(dset,'waveB2Value',nstart_2d,ncount_2d);
%   mask_2d = repmat(mask_1d,1,size(freq,1))';
%   if strcmp(MopQCflag(1:2),'on')
%    na0(mask_2d~=1) = NaN;
%   end
%   % na1(mask_2d~=1) = NaN;
%   % nb1(mask_2d~=1) = NaN;
%   % na2(mask_2d~=1) = NaN;
%   % nb2(mask_2d~=1) = NaN;
% 
% end
% 
% % ecmwf forecast is in 3 hourly steps. Interpolate to hourly
% nptime=(nptime(1):1/24:nptime(end))';
% na0=interp1(1:size(na0,2),na0',1:1/3:size(na0,2))';
% 
% % forecast time series starts before the current time
% % find overlap between nowcast and forecast data
% [owavetime,iht,ift]=intersect(moptime,nptime); % overlaping hindcast and forecast times
% 
% % make energy bias correction to forecast
% waveE=a0'*bw;
% fwaveE=na0'*bw;
% Ebias=mean(waveE(iht)./fwaveE(ift)); % calculate Hs bias
% % % make wave height bias correction
% na0=na0*Ebias;
% 
% % remove overlap time from forecast
% %wtime2(ift)=[];
% na0(:,ift)=[];
% nptime(ift)=[];
% 
% % combine nowcast and forecast matrices
% moptime=[moptime' nptime']';
% a0=[a0 na0];
% fidx=(numel(moptime)-numel(nptime)):numel(moptime);
% % a1=[a1 na1];
% % b1=[b1 nb1];
% % a2=[a2 na2];
% % b2=[b2 nb2];
% 
% return
% end
% %% ---------------------------------------------------------------------
% 
% function [L,C,Cg]=LinearDispersion(frequency,depth)
% 
% %% Linear Dispersion Equation Solver 
% 
% % Input: wave frequency (Hz), and water depth (m)
% 
% % Returns: wavelength L (m) , phase speed c (m/s), 
% %          and group speed Cg (m/s) of a linear wave.
% 
% % Solves the linear dispersion relationship for surface gravity waves
% 
% %    (2pi/frequency)^2=g*(2pi/L)*tanh(depth*2pi/L)
% 
% %  for L, where g= gravitational acceleration (9.81m/s^2), 
% %   using an algorithm derived by C.S. Wu.
% 
% % Wu, Chung-Shang & Thornton, Ed. (1986). Wave Numbers of Linear Progressive 
% % Waves. J WATERW PORT COAST OC-ASCE. 112. 10.1061/(ASCE)0733-950X(1986)112:4(536). 
% 
% % This is one of many published approximate solutions. 
% 
% % For a detailed derivation and explanation of the linear dispersion
% %  relationship see: http://falk.ucsd.edu/pdf/WavesLecture02_211A.pdf 
% 
% %% Apply the C.S. Wu algorithm to get radian wavenumber
% 
% a=4.0243*depth*frequency.^2;
% yhat=a.*(1+1.26*exp((-1.84).*a));t=exp((-2).*yhat);
% aka=a;aka(a >= 1)=a(a >= 1).*(1+2*t(a >= 1).*(1+t(a >= 1)));
% aka(a < 1)=sqrt(a(a < 1)).*(1+(a(a < 1)./6).*(1+a(a < 1)./5));
% k=abs(aka./depth); % radian wavenumber 2pi/L
% 
% % linear wave parameters
% 
% L=2*pi./k; % wavelength in meters
% C=L.*frequency; % wave speed in m/s
% % group velocity in m/s
% sigma=2*depth.*k;
% Cg=pi*(frequency./k).*(1+sigma./sinh(sigma));
% 
% end