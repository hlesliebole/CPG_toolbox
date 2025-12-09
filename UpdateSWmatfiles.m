%  updates existing M*SW.mat file

%  make survey waves SW struct array mat file based
%   on survey dates in the survey averaged data SA mat file
clear all
%MopNumber=503;
for MopNumber=616:626 %495:626

fprintf('%i of 626\n',MopNumber)
stn = ['D' num2str(MopNumber,'%4.4i')];
SAmatfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(SAmatfile)

% load exist wave mat file
SWmatfile=['M' num2str(MopNumber,'%5.5i') 'SW.mat'];
load(SWmatfile)
% retain as old SW data
OW=SW;
clear SW

% sort SA by date
T=struct2table(SA);
sortedT = sortrows(T, 'Datenum');
SA=table2struct(sortedT)';
clear T sortedT;

urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';  
% Set 'base' THREDDS URL and pieces to concatenate with user-defined station number (set above)
urlend = '_hindcast.nc'; % change this to _forecast.nc to make forecast prediction
dsurl=strcat(urlbase,stn,urlend);

% read in netcdf wave hindcast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime');

hs=ncread(dsurl,'waveHs');
sxy=ncread(dsurl,'waveSxy');
wtime=wavetime;
whs=hs;
wsxy=sxy;

urlend = '_nowcast.nc'; % change this to _forecast.nc to make forecast prediction
dsurl=strcat(urlbase,stn,urlend);

% read in netcdf wave hindcast times 
wavetime=ncread(dsurl,'waveTime');
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime');

wavetime=vertcat(wtime,wavetime);

hs=ncread(dsurl,'waveHs');
sxy=ncread(dsurl,'waveSxy');
wavehs=vertcat(whs,hs);
waveSxy=vertcat(wsxy,sxy);
% wave energies
E=(wavehs/4).^2;

% no data for first date
SW(1).Mopnum=MopNumber;
SW(1).Datenum=SA(1).Datenum;
SW(1).NetHours=0;
SW(1).Emean=0; 
SW(1).Emedian=0;
SW(1).Emax=0;
SW(1).E90q=0;
SW(1).E95q=0;
SW(1).E98q=0;
SW(1).SxyMean=0;
SW(1).SxyMedian=0;
SW(1).SxyMin=0;
SW(1).SxyMax=0;
SW(1).SxyNet=0;

% loop through remaining survey dates and calculate wave stats
%  since previous survey
for ns=2:size(SA,2)
    
% min,max,mean between surveys
SW(ns).Mopnum=MopNumber;
SW(ns).Datenum=SA(ns).Datenum;

% see if old data exist for this survey time AND date of previous
%  survey matches as well.
% check for old match to survey date
old=find([OW.Datenum] == SA(ns).Datenum);
% check if previous survey date matches as well
if ~isempty(old) 
    if(old(1) > 1 && OW(old(1)-1).Datenum ~= SA(ns-1).Datenum)
        old=[]; % no match, reset old to empty
    end
end

if ~isempty(old) % copy over old data 
        SW(ns).NetHours=OW(old).NetHours;
        SW(ns).Emean=OW(old).Emean;
        SW(ns).Emedian=OW(old).Emedian;
        SW(ns).Emax=OW(old).Emax;
        SW(ns).E90q=OW(old).E90q;
        SW(ns).E95q=OW(old).E95q;
        SW(ns).E98q=OW(old).E98q;
        SW(ns).SxyMean=OW(old).SxyMean;
        SW(ns).SxyMedian=OW(old).SxyMedian;
        SW(ns).SxyMin=OW(old).SxyMin;
        SW(ns).SxyMax=OW(old).SxyMax;
        SW(ns).SxyNet=OW(old).SxyNet;
else
  % calculate new data
  fprintf('adding %s \n',datestr(SW(ns).Datenum))
i=find(wavetime >= datestr(SA(ns-1).Datenum) & wavetime <= datestr(SA(ns).Datenum+1));
%fprintf('%i %i %s to %s\n',ns,length(i),datestr(SA(ns-1).Datenum),datestr(SA(ns).Datenum));
if ~isempty(i)
    SW(ns).NetHours=length(i);

% min,max,mean between surveys
SW(ns).Emean=nanmean(E(i)); 
SW(ns).Emedian=nanmedian(E(i));
SW(ns).Emax=max(E(i));

Eq=quantile(E(i),[.90 .95 .98]);
SW(ns).E90q=Eq(1);
SW(ns).E95q=Eq(2);
SW(ns).E98q=Eq(3);

SW(ns).SxyMean=nanmean(waveSxy(i));
SW(ns).SxyMedian=nanmedian(waveSxy(i));
SW(ns).SxyMin=nanmin(waveSxy(i));
SW(ns).SxyMax=nanmax(waveSxy(i));
SW(ns).SxyNet=nansum(waveSxy(i));
end

end

end

SWmatfile=['M' num2str(MopNumber,'%5.5i') 'SW.mat'];
eval(['save ' SWmatfile ' SW']);

end
