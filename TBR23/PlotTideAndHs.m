%figure('position',[63   306   810   466])
Mopnum=582;

dtm=datenum(2023,4,5,1,0,0);
dtm2=datenum(2024,1,31,1,0,0);

stn=['D0' num2str(Mopnum)];
% create url pathname to CDIP thredds netcdf nowcast file
urlbase=...
 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';
urlend='_ecmwf_fc.nc';% get the MOP nowcast file
dsurl = strcat(urlbase,stn,urlend);% put the url string together
fwavetime=ncread(dsurl,'waveTime');% read in netcdf wave times 
% convert mop unix time to local time
fwavetime=datetime(fwavetime,'ConvertFrom','posixTime',...
    'TimeZone','America/Los_Angeles');
% get energy densities and f bands
% a0=ncread(dsurl,'waveEnergyDensity');
% f=double(ncread(dsurl,'waveFrequency'));
% bw=ncread(dsurl,'waveBandwidth');
fHs=ncread(dsurl,'waveHs');
fTp=ncread(dsurl,'waveTp');

stn=['D0' num2str(Mopnum)];
% create url pathname to CDIP thredds netcdf nowcast file
urlbase=...
 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';
urlend='_nowcast.nc';% get the MOP nowcast file
dsurl = strcat(urlbase,stn,urlend);% put the url string together
wavetime=ncread(dsurl,'waveTime');% read in netcdf wave times 
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime',...
    'TimeZone','America/Los_Angeles');
% get energy densities and f bands
% a0=ncread(dsurl,'waveEnergyDensity');
% f=double(ncread(dsurl,'waveFrequency'));
% bw=ncread(dsurl,'waveBandwidth');
Hs=ncread(dsurl,'waveHs');
Tp=ncread(dsurl,'waveTp');

% convert closest local noon mop wave record into energies in m^2
[~,wstart]=min(abs(dtm-datenum(wavetime)));wstart=wstart+12;
[~,wend]=min(abs(dtm2-datenum(wavetime)));wend=wend+12;
idx=wstart:wend;

% [tidetime,tidehgt]=...
%     getztide2(datenum(wavetime(idx(1))),datenum(wavetime(idx(end))),...
%     'gmt','navd','predictions');
[tidetime,tidehgt]=...
    getztide2(dtm,dtm2,...
    'gmt','navd','predictions');
tidetime=datetime(datenum(tidetime),'convertfrom','datenum','TimeZone','America/Los_Angeles');

%ndx=find(ismember(wavetime,tidetime));
ndx=idx;
%twl=tidehgt+Hs(ndx);

%a2=axes('position',[.05 .05 .9 .25]);
pt=plot(tidetime,tidehgt,'-','color',[.7 .7 .7],'linewidth',2,'DisplayName','Tide');hold on; %navd88 tide elev

pf=plot(fwavetime,fHs,'-','color',[0 .7 0],'DisplayName','Forecast Hs','linewidth',2);hold on;set(gca,'fontsize',12);
plot(wavetime(wstart:end),Hs(wstart:end),'k-','DisplayName','Hs');hold on;set(gca,'fontsize',12);
p1=plot(wavetime(wstart:end),Hs(wstart:end),'k-','linewidth',2,'DisplayName','Hs');
%p3=plot(tidetime,twl,'m-','linewidth',2,'DisplayName','Tide+Hs');
set(gca,'xlim',[wavetime(wstart) wavetime(end)]);
ylabel('NAVD88 Tide Elev and Hs (m)');set(gca,'ylim',[-.7 2.2])
yyaxis right
pt=plot(fwavetime,fTp,'c.','markersize',10,'DisplayName','Forecast Tp');
p2=plot(wavetime(wstart:end),Tp(wstart:end),'r.','markersize',10,'DisplayName','Tp');

grid on;%set(a1,'xdir','reverse','color',[.8 .8 .8]);
ylabel('Tp (s)');set(gca,'ycolor','r')
xlabel('Date');title(['Mop ' num2str(Mopnum) ' Tides and Waves']);
%legend([pt p1 p2],'location','northwest')
set(gca,'xlim',[tidetime(1) tidetime(end)])

%makepng('TidesAndWaves.png')