function [freq,bw,depth,moptime,fidx,a0] = GetMOPfreqSpectra(MopName,dstart,dend,MopQCflag)

% returns spectral data from the nowcast and ecmwf forecast thredds
%  netcdf files for the mop specified by MopID and dstart to dend datenumber time range

% create urls for MopID on thredds server
urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';  
dset=strcat(urlbase,MopName,'_nowcast.nc'); % nowcast url
dset2=strcat(urlbase,MopName,'_ecmwf_fc.nc'); % forecast url
% Display variable names in the netcdf file in the command window
%ncdisp(dset,'/','min')

% load nowcast time variable, find time indices needed
  wtime = ncread(dset,'waveTime');
  % convert posix time to datenumber
  wtime = double(wtime)/86400 + datenum(1970,1,1);
  datestr(wtime(end))

  % load forecast time variable, find time indices needed
  wtime2 = ncread(dset2,'waveTime');
  % convert posix time to datenumber
  wtime2 = double(wtime2)/86400 + datenum(1970,1,1);

%  option for datetime variable
% wavedatetime=datetime(wtime,'ConvertFrom','posixTime');

%% get spectral data from the nowcast file

a0=[];a1=[];b1=[];a2=[];b2=[];moptime=[];
% check if at least some of the data is in the nowcast file
if dstart <= wtime(end) 

  dindices = find(wtime >= dstart & wtime <= dend);
  dlength = size(dindices,1);
  nstart_1d = dindices(1);
  ncount_1d = dlength;
  moptime = wtime(dindices(1):dindices(ncount_1d));
  % load qc flags; records good where qcflag = 1 
  qcflag = ncread(dset,'waveFlagPrimary',nstart_1d,ncount_1d);
  mask_1d = qcflag==1;
  
   % load 2d var (a1), mask off before plotting
  freq = ncread(dset,'waveFrequency');
  bw = ncread(dset,'waveBandwidth');
  depth=ncread(dset,'metaWaterDepth');
  nstart_2d = [1, dindices(1)];
  ncount_2d = [size(freq,1), dlength];
  a0 = ncread(dset,'waveEnergyDensity',nstart_2d,ncount_2d);
  % a1 = ncread(dset,'waveA1Value',nstart_2d,ncount_2d);
  % a2 = ncread(dset,'waveA2Value',nstart_2d,ncount_2d);
  % b1 = ncread(dset,'waveB1Value',nstart_2d,ncount_2d);
  % b2 = ncread(dset,'waveB2Value',nstart_2d,ncount_2d);
  mask_2d = repmat(mask_1d,1,size(freq,1))';
  if strcmp(MopQCflag(1:2),'on')
   a0(mask_2d~=1) = NaN;
  end
  % a1(mask_2d~=1) = NaN;
  % b1(mask_2d~=1) = NaN;
  % a2(mask_2d~=1) = NaN;
  % b2(mask_2d~=1) = NaN;

end

%% get forecast file data

na0=[];na1=[];nb1=[];na2=[];nb2=[];nptime=[];

if dend >= wtime2(1)

  dindices = find(wtime2 >= dstart & wtime2 <= dend);
  dlength = size(dindices,1);
  nstart_1d = dindices(1);
  ncount_1d = dlength;
  nptime = wtime2(dindices(1):dindices(ncount_1d));

  % load qc flags; records good where qcflag = 1 
  qcflag = ncread(dset,'waveFlagPrimary',nstart_1d,ncount_1d);
  mask_1d = qcflag==1;
  
   % load 2d var (a1), mask off before plotting
  freq = ncread(dset,'waveFrequency');
  bw = ncread(dset,'waveBandwidth');
  depth=ncread(dset,'metaWaterDepth');
  nstart_2d = [1, dindices(1)];
  ncount_2d = [size(freq,1), dlength];
  na0 = ncread(dset,'waveEnergyDensity',nstart_2d,ncount_2d);
  % na1 = ncread(dset,'waveA1Value',nstart_2d,ncount_2d);
  % na2 = ncread(dset,'waveA2Value',nstart_2d,ncount_2d);
  % nb1 = ncread(dset,'waveB1Value',nstart_2d,ncount_2d);
  % nb2 = ncread(dset,'waveB2Value',nstart_2d,ncount_2d);
  mask_2d = repmat(mask_1d,1,size(freq,1))';
  if strcmp(MopQCflag(1:2),'on')
   na0(mask_2d~=1) = NaN;
  end
  % na1(mask_2d~=1) = NaN;
  % nb1(mask_2d~=1) = NaN;
  % na2(mask_2d~=1) = NaN;
  % nb2(mask_2d~=1) = NaN;

end

% ecmwf forecast is in 3 hourly steps. Interpolate to hourly
nptime=(nptime(1):1/24:nptime(end))';
na0=interp1(1:size(na0,2),na0',1:1/3:size(na0,2))';

% forecast time series starts before the current time
% find overlap between nowcast and forecast data
[owavetime,iht,ift]=intersect(moptime,nptime); % overlaping hindcast and forecast times

% make energy bias correction to forecast
waveE=a0'*bw;
fwaveE=na0'*bw;
Ebias=mean(waveE(iht)./fwaveE(ift)); % calculate Hs bias
% % make wave height bias correction
na0=na0*Ebias;

% remove overlap time from forecast
%wtime2(ift)=[];
na0(:,ift)=[];
nptime(ift)=[];

% combine nowcast and forecast matrices
moptime=[moptime' nptime']';
a0=[a0 na0];
fidx=(numel(moptime)-numel(nptime)):numel(moptime);
% a1=[a1 na1];
% b1=[b1 nb1];
% a2=[a2 na2];
% b2=[b2 nb2];

return
end