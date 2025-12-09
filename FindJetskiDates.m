% Example code to make a plot of the survey date vs the 
%  time between jetski surveys

Mopnum=667; % center of Cardiff
%Mopnum=584; % Torrey Pines Parking Lot

% Paths to CPGMOP data files and the MOPS/toolbox
mpath='/volumes/group/MOPS/'; % reefbreak on a mac

% load the mat file containing processed survey profiles
eval(['load ' mpath 'M00' num2str(Mopnum,'%3.3i') 'SM.mat']);

% find surveys with original data file name that included the word jumbo
jumbo=find(contains({SM.File},'umbo'));

% Some files with the jumbo name do not include jetski data
%  so find jumbo survey subset with profile data below -3m navd88
m=0;
jetski=[];
for j=1:length(jumbo)
    fprintf('%s %5.1f\n',datestr(SM(jumbo(j)).Datenum),min(SM(jumbo(j)).Z1Dmean));
    if( min(SM(jumbo(j)).Z1Dmean) < -3 )
        m=m+1;
        jetski(m)=jumbo(j);
    end       
end

figure('position',[66         186        1149         553]);
x=datetime([SM(jetski).Datenum],'convertfrom','datenum');
x(1)=[]; % drop first date for plotting days between surveys
plot(x,diff([SM(jetski).Datenum]),'.-','markersize',20,'linewidth',2);
hold on;
set(gca,'yscale','log','ylim',[1 500],'fontsize',14,'linewidth',2);
grid on
set(gca,'xlim',[datetime(year(x(1))-1,1,0,0,0,0) datetime(year(x(end))+1,1,0,0,0,0)]);
xl=get(gca,'xlim');
text(xl(2),1,' - 1 Day','verticalalign','middle','fontsize',14)
text(xl(2),4,' - 4 Days','verticalalign','middle','fontsize',14)
text(xl(2),7,' - 1 week','verticalalign','middle','fontsize',14)
text(xl(2),14,' - 2 weeks','verticalalign','middle','fontsize',14)
text(xl(2),30,' - 1 month','verticalalign','middle','fontsize',14)
text(xl(2),91,' - 3 months','verticalalign','middle','fontsize',14)
text(xl(2),182,' - 6 months','verticalalign','middle','fontsize',14)
xlabel('Date')
ylabel('Days from Previous Jetski Survey')
title([' Mop ' num2str(Mopnum) ' Historical Jetski Surveys'],'fontsize',18)
set(gca,'xtick',datetime(datenum(year(x(1))-1:year(x(end))+1,1,1),'convertfrom','datenum')');

% make png image usinfg mopnum in name
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','-r300','-loose',['Mop' num2str(Mopnum,'%3.3i') 'JetskiDates.png']);

fprintf('\n%s\n%s\n','Plot written to:',['Mop' num2str(Mopnum,'%3.3i') 'JetskiDates.png']);
