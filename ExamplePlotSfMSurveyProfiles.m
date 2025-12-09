% Example code to plot all the SfM drone survey
%  profiles at a specific Mop

addpath /volumes/group/Mops
addpath /volumes/group/Mops/toolbox

MopNumber=582; % select Mop

load(['M00' num2str(MopNumber,'%3.3i') 'SM.mat']) % load Mop SM file
ndx=find(strcmp({SM.Source}','SfMdrone')); % find SfM drone surveys

% plot symbols to use randomly 
stype=cell(numel(ndx),1);
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

figure('position',[1          55        1412         742]);
m=0;
for n=ndx'   
     m=m+1;
     t(m)=SM(n).Datenum; % date vector
     % plot the transect profile with randon color and symbol
         p(m)=plot(SM(n).X1D,SM(n).Z1Dtransect,'LineStyle','-',...
        'Marker',all_marks{mod(1+round(11*rand),13)},...
        'DisplayName',[datestr(SM(n).Datenum) ' ' SM(n).Source],...
        'linewidth',2);hold on;        
end

grid on;
set(gca,'fontsize',12);xlabel('Distance From Mop Back Beach Point (m)');ylabel('Elevation (m, NAVD88)')
set(gca,'ylim',[-1 3.5],'xdir','reverse','fontsize',14);

title([{['Mop ' num2str(MopNumber)]},{'SfM Drone Surveys'}],...
   'fontsize',16);
xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);  

legend(p,'location','eastoutside','fontsize',12,'numcolumns',ceil(numel(ndx)/60))
makepng(['SfMmop' num2str(MopNumber) 'Profiles.png'])
