% MatchMyBeachFace.m

% This script finds the jumbo full mean profile at a Mop
% that best matches a specified subaerial survey's date 
% beach face (eg. 0 to 2m navd88).

% Idea is to make a best guess at the full profile based
%  on, eg. the lastest subaerial beach profile.


% settings
% --------

% link to cpg mop system, edit as needed
addpath /volumes/group/MOPS
addpath /volumes/group/MOPS/toolbox

% Mop Number
MopNumber=582;
% Beachface elevation (navd88) range to fit jumbos to 
zbfmin=0.5;
zbfmax=1.0;
% Subaerial Survey Date Number to match jumbos to.
% Set = 0 if you want it to fit jumbos to the most
%  recent not-jumbo survey
SSdate=0; %SSdate=datenum(2021,10,11); % 10/11/2021 prestorm example

% load the CPG MOP survey morpho (SM) struct array
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');
% if isempty(SM(end).File)
%     SM(end).File='/volumes/group/ig8_wheel/20211013_00581_00583_Torrey_wheel/filtered_clean20211013.llnezts.navd88';
%     SM(end).Datenum=datenum(2021,10,13);
% end
% if isempty(SM(end-1).File)
%     SM(end-1).File='/volumes/group/ig8_wheel/20211012_00581_00583_Torrey_wheel/filtered_clean20211012.llnezts.navd88';
%     SM(end-1).Datenum=datenum(2021,10,12);
% end
% if isempty(SM(end-2).File)
%     SM(end-2).File='/volumes/group/ig8_wheel/20211011_00581_00583_Torrey_wheel/filtered_clean20211011.llnezts.navd88';
%     SM(end-2).Datenum=datenum(2021,10,11);
% end
% save survey morpho results in SM file for Mop number
%save(['M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM');

% find the indexes of all the jumbos
JumboIndexes=find(contains({SM.File}, 'jumbo','IgnoreCase',true)==1); 

% find the index of the target subaerial survey date to match
if SSdate > 0
 TargetIndex=find([SM.Datenum] == SSdate);
else
  % find the index of the most recent survey that is NOT a jumbo
  TargetIndex=find(contains({SM.File}, 'jumbo','IgnoreCase',true)==0,1,'last'); 
end

% loop through jumbos and get rmse fit to the beach face over the
%  defined navd88 elevation range

% xshore range of target beach face elevations

xrange=find(SM(TargetIndex).Z1Dmean >= zbfmin &...
    SM(TargetIndex).Z1Dmean <= zbfmax);
xrange=min(xrange):max(xrange); % make sure the xshore range is contiguous

m=0;
rmse=[];ngaps=[];
for n=JumboIndexes
    %jxrange=find(SM(n).X1D == xrange(1)):find(SM(n).X1D == xrange(end));
    jxrange=xrange;
    err=SM(n).Z1Dmean(jxrange)-SM(TargetIndex).Z1Dmean(xrange);
    m=m+1;
    rmse(m)=sqrt(nanmean(err.^2));
    ngaps(m)=sum(isnan(err));
end

[rmin,imin]=min(rmse);

figure;
p(1)=plot(SM(JumboIndexes(imin)).X1D,SM(JumboIndexes(imin)).Z1Dmean,'b-',...
    'DisplayName',datestr(SM(JumboIndexes(imin)).Datenum,'mm/dd/yyyy')); 
hold on;
p(2)=plot(SM(TargetIndex).X1D,SM(TargetIndex).Z1Dmean,'r-',...
    'DisplayName',datestr(SM(TargetIndex).Datenum,'mm/dd/yyyy'));
p(3)=plot(SM(JumboIndexes(end)).X1D,SM(JumboIndexes(end)).Z1Dmean,'g-',...
    'DisplayName',datestr(SM(JumboIndexes(end)).Datenum,'mm/dd/yyyy'));
set(gca,'xdir','reverse');
grid on;
xl=get(gca,'xlim');yl=get(gca,'ylim');
xlabel('Xshore Distance (m)');
ylabel('Date');
zlabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(SM(1).Mopnum) ' Area Mean Xshore Profiles']);
set(gca,'fontsize',14);
legend(p,'location','northwest');
makepng('BestFitJumboToLastBeachface.png')
JumboIndexes(imin)









