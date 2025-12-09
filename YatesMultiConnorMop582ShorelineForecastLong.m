% Example script to predict the change in the MSL shoreline location
%  since the last survey.  Uses Yates et al equations. Uses Thredds
%  MOPS nowcast and forecast netcdf files.

% makes bias-corrected ww3-driven Hs forecast based on overlap between 
% (buoy-driven) MOP hindcast/nowcast and (ww3-driven) MOP forecast

% The msl location is plotted in meters from a MOP's back beach location.
clear all
%close all

%-------------------------------
%   Script Settings
%--------------------------------
%close all
% MOP Setting 
%-----------------------------------------------------
% plot shoreline location and wave height time series 
%------------------------------------------------------
%figure('position',[166   120   962   654]);
figure('position',[20  20   962*1.45   654*1.2]);

col=flipud(colormap(jet(8)));col(4,:)=0.8*col(4,:);

% Initial MSL Shoreline Settings

% 1. Estimated mean annual MSL distance from the back beach point 

% This is the global GM.X1D xshore distance of MSL from the GM.Z1Dmean
%  profile in M00582GM.mat
meanS=43; 

% 2. Most recent processed survey date [ SM(end).Datenum ] in M00582SM.mat 
%survdate=datetime(2021,10,20,18,0,0,'TimeZone','America/Los_Angeles');
survdate=datetime(2021,10,28,18,0,0,'TimeZone','America/Los_Angeles');

% 3. Estimated distance to msl shoreline in most recent survey

% SM(end).X1D xshore distance to MSL (0.774 navd88) in SM(end).Z1Dmean
%survS=38.2; 
survS=31.65; 

% Yates MSL Model Settings
%  using model coefficients from torrey pines
b=0.07;
a=-0.0045;
Cminus=-1.38;
Cplus=-1.16;
%Cplus=0; % turn accretion off

%--------------------------------------------------------------------------
%  ------ Create shoreline hindcast/nowcast since the last survey ---------
%--------------------------------------------------------------------------
for MopNumber=578:584
stn=['D0' num2str(MopNumber)]; 

% create url pathname to CDIP thredds netcdf directory
urlbase=...
 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';
urlend='_nowcast.nc';% get the MOP nowcast file
dsurl = strcat(urlbase,stn,urlend);% put the url string together
wavehs=ncread(dsurl,'waveHs'); % read in netcdf wave hs nowcast
E=(wavehs/4).^2;% convert to energy
wavetime=ncread(dsurl,'waveTime');% read in netcdf wave times 
% convert mop unix time to local time
wavetime=datetime(wavetime,'ConvertFrom','posixTime',...
    'TimeZone','America/Los_Angeles');


%  ------ run yates model from most recent survey date to the latest
%         mop nowcast prediction hour  ---------------------------


t1=survdate;t2=wavetime(end);  % start and end times

% preallocate array large enough to hold results
S=zeros(1,length(E));

% set the first shoreline location to the latest survey location
%  minus the annual mean shoreline location
S(find(wavetime == t1)-1)=survS-meanS;

% step through hourly wave data and change the shoreline location

for i=find(wavetime == t1):find(wavetime == t2)
    Eeq=a*S(i-1)+b; % Yates eq. 4 equilibrium E
    deltaE=E(i)-Eeq; % Yates eq. 3 current diff from equilibrium E
    if( deltaE > 0)   
        dSdt=Cminus*sqrt(E(i))*deltaE; % Yates eq. 2 erosion
    else     
        dSdt=Cplus*sqrt(E(i))*deltaE; % Yates eq. 2 accretion
    end
    
    S(i)=S(i-1)+dSdt; % change location
end


% isolate mop hindcast data time range to plot since last survey
i=find(wavetime >= t1 & wavetime <= t2);


%--------------------------------------------------------------------------
%  ------ Create shoreline forecast based on mop wave forecast ---------
%--------------------------------------------------------------------------

urlend='_forecast.nc';% get the MOP forecast file
dsurl = strcat(urlbase,stn,urlend);% put the url string together
fwavehs=ncread(dsurl,'waveHs'); % read in netcdf wave hs nowcast
fE=(fwavehs/4).^2;% convert to energy
fwavetime=ncread(dsurl,'waveTime');% read in netcdf wave times 
% convert mop unix time to local time
fwavetime=datetime(fwavetime,'ConvertFrom','posixTime',...
    'TimeZone','America/Los_Angeles');

% make bias corrected Hs based on overlap between (buoy-driven) MOP 
% hindcast/nowcast and (ww3-driven) MOP forecast

[owavetime,iht,ift]=intersect(wavetime,fwavetime); % overlaping hindcast and forecast times
Hsbias=mean(wavehs(iht)./fwavehs(ift));
cfwavehs=fwavehs*Hsbias;cfE=(cfwavehs/4).^2;% convert to corrected energy

% make Yates shoreline change forecast from the last location in the 
%  hindcast/forecast shoreline location time series

t1=owavetime(end);t2=fwavetime(end);  % start and end times

% preallocate array large enough to hold results
fS=zeros(1,length(fE));

% set the first shoreline location to the last hindcast/nowcast location
%  minus the annual mean shoreline location
fS(find(fwavetime == t1)-1)=S(iht(end));

% Note: mop wave forecasts have a 3 hour time step, so need to multiply
%  the resulting Yates shoreline change amounts (Cplus and Cminus based
%  on hourly data) by a factor of 3

% step through 3-hourly hourly corrected forecast wave data and change 
%  the shoreline location

for j=find(fwavetime == t1):find(fwavetime == t2)
    Eeq=a*fS(j-1)+b; % Yates eq. 4 equilibrium E
    deltaE=cfE(j)-Eeq; % Yates eq. 3 current diff from equilibrium E
    if( deltaE > 0)   
        dSdt=Cminus*sqrt(cfE(j))*deltaE; % Yates eq. 2 erosion
    else     
        dSdt=Cplus*sqrt(cfE(j))*deltaE; % Yates eq. 2 accretion
    end
    
    fS(j)=fS(j-1)+dSdt*3; % change location * 3hourly factor of 3
end

% plot wave forecast
subplot(2,1,1);
if MopNumber == 582
ph(MopNumber-577)=plot(wavetime(i),wavehs(i),'k-','linewidth',2);hold on;
set(ph(MopNumber-577),'displayname',[num2str(MopNumber)]); 
else
ph(MopNumber-577)=plot(wavetime(i),wavehs(i),'-','linewidth',1);hold on;
set(ph(MopNumber-577),'color',col(MopNumber-577,:),'displayname',[num2str(MopNumber)]); 
end
%col(MopNumber-577,:)=get(ph(MopNumber-577),'color');


if MopNumber == 584
pf1=plot(fwavetime,fwavehs,'r-','linewidth',1,'displayname','Hs Forecast (ww3-driven)');%datetick('x','mmm dd');
else
pf=plot(fwavetime,fwavehs,'r-','linewidth',1);
end
grid on;hold on;%xlim([fwavetime(1) fwavetime(end)]);
%xlim([survdate fwavetime(end)]);
if MopNumber == 584
pcf1=plot(fwavetime,cfwavehs,'r--','linewidth',1,'displayname','Bias-Corrected Hs Forecast');
else
pcf=plot(fwavetime,cfwavehs,'r--','linewidth',1);
end
yl=get(gca,'ylim');set(gca,'ylim',[0 2.5]);%set(gca,'ylim',[yl(1) yl(2)+.25*(diff(yl))]);
%legend([ph pf pcf],'location','northwest');
%legend('location','northwest');
ylabel('Forecast Hs(m)');xlabel('Local Time');
title(['MOPs 578-584 Buoy-drive Wave Hindcast and WW3-driven Wave Forecast']);
set(gca,'xlim',[dateshift(survdate,'start','day') dateshift(fwavetime(end),'end','day')])
xl=get(gca,'xlim');set(gca,'Xtick',xl(1):xl(2))
datetick('x','mm/dd','keepticks')
if MopNumber == 584
%legend([ph pf1 pcf1],'location','northeast');
legend(ph,'location','northwest');
end

% plot msl forecast
subplot(2,1,2);

yyaxis left
[tidetime,tidehgt]=getztide2(datenum(wavetime(1)),datenum(fwavetime(end)),'lst','msl','predictions');
xtide=datenum(tidetime);
plot(datetime(xtide,'convertfrom','datenum','TimeZone','America/Los_Angeles'),tidehgt,...
    '-','color',[.8 .8 .8],'linewidth',2)
ylabel('Tide (m, MSL)','color','k');set(gca,'ycolor','k');


yyaxis right


if MopNumber == 582
    phf=plot(wavetime(i),S(i)+meanS,'k-','linewidth',2);hold on;
else
phf=plot(wavetime(i),S(i)+meanS,'-','linewidth',2,'color',col(MopNumber-577,:));hold on;
end
ii=find(fwavetime >= t1 & fwavetime <= t2);
phf=plot(fwavetime(ii),fS(ii)+meanS,'r-','linewidth',2);hold on;

grid on;hold on;%xlim([fwavetime(1) fwavetime(end)]);
xlim([dateshift(survdate,'start','day') dateshift(fwavetime(end),'end','day')]);
pm=plot([survdate fwavetime(end)],[meanS meanS],'k--','linewidth',2); 
%yl=get(gca,'ylim');set(gca,'ylim',[yl(1) yl(2)+.2*(diff(yl))]);

title({['Hindcast and Forecast MSL Shoreline Change since ' datestr(survdate) ' Survey'];...
    'All MOP xMSL distances adjusted relative to the MOP 582 starting xMSL value of 38.2m' });
ylabel('Meters from Back Beach Point');xlabel('Local Time');
% adjust y axis limits to show mean shoreline if necessary
% yl=get(gca,'ylim');if(yl(1) == meanS);yl(1)=meanS-5;end
% if(yl(2) == meanS);yl(2)=meanS+5;end
% set(gca,'ylim',yl);
% match date range of wave forecast plot

set(gca,'xlim',[dateshift(survdate,'start','day') dateshift(fwavetime(end),'end','day')])
xl=get(gca,'xlim');set(gca,'Xtick',xl(1):xl(2))
datetick('x','mm/dd','keepticks')

end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------


% get recent 582 xMSL locations from CPG MOP database
addpath /volumes/group/mops

for MopNumber=578:584
%MopNumber=582;

% load morpho data for the mop
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');

% find survey dates starting from forecast survdate beginning
idx=find([SM.Datenum] >= datenum(survdate-1));

%nn=0;
iflag=0;
for n=1:numel(idx) 
    m=idx(n);
%     if n==1;SM(m).Datenum=datenum(2021,10,11);end
%     if n==2;SM(m).Datenum=datenum(2021,10,12);end
%     if n==3;SM(m).Datenum=datenum(2021,10,13);end
    if ~isnan(min(SM(m).Z1Dmean))
    %nn=nn+1;
     z=SM(m).Z1Dmean;z(z <0)=NaN;
     xMSL=intersections([SM(m).X1D(1) SM(m).X1D(end)],[0.774 0.774],SM(m).X1D,SM(m).Z1Dmean);
     if ~isempty(xMSL);xMSL=xMSL(end);end
     otime=datetime(SM(m).Datenum+.75,'convertfrom','datenum','TimeZone','America/Los_Angeles');
     
     if iflag == 0 && ~isempty(xMSL); dx=survS-xMSL; iflag=1;end
     if ~isempty(xMSL)
         xMSL=xMSL(end);
     if MopNumber == 582
         pobs(MopNumber-577)=plot(otime,xMSL+dx,'k.','markersize',20,'displayname',[num2str(MopNumber)]);
     else
     pobs(MopNumber-577)=plot(otime,xMSL+dx,'o','markersize',5,'linewidth',2,'color',col(MopNumber-577,:),'displayname',[num2str(MopNumber)]);
     end
     end
 
%     p(nn)=plot(SM(m).X1D,z,'-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yyyy') ' xMSL=' num2str(xMSL,'%4.1f') 'm']);hold on;
    end
end

end
% legend([pm phf pf pobs],'Mean Shoreline Location','MSL Hindcast/Nowcast',...
%     'Bias-Corrected MSL Forecast','Connor GPS Wheel of Fortune','location','northwest');
legend(pobs,'location','northwest');
set(gca,'ylim',[20 50])

% make png file

pfile=['MopMultishorelineForecast.png'];
print(gcf,'-dpng','-r300','-loose',pfile);

