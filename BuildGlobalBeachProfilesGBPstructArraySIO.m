%%
%  Code to read CPG Mops *SA.mat files and create a Global Subaerial Beach Profile
%   struct array GBP
clear all

addpath /Users/William/Desktop/MOPS
addpath /Users/William/Desktop/MOPS/toolbox

nmop=0;
MopStart=507; % Scripps Pier ; %520 Scripps Canyon
MopEnd=514; % TP north end; 683; % San Elijo Inlet

% id mops that require additional qc before 
%  making mean profiles
mopedit=[];%[534 535 555 556];

for MopNumber=MopStart:MopEnd 
    fprintf('%i\n',MopNumber)
    nmop=nmop+1;

matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');

% make non-gridded transects
for n=1:size(SA,2)
    [x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
    SA(n).X1D=x1d;
    SA(n).Z1Dtransect=z1di;
end

% identify jumbos
jumbo=find(contains({SA.File},'umbo'));

% % reduce local SA struct array to just jumbo surveys
%SA=SA(jumbo);

% % edit SA data for mops with known qc issues
% if ismember(MopNumber,mopedit)
%     SA=EditSA(SA);
% end

%  Get the beach year for each jumbo by adding 92 days to the 
%  surveys datenum to push Oct 1 into the next year
BeachYr=year(datetime([SA.Datenum]+92,'convertfrom','datenum'));

% Get the months of the surveys
SurvMon=month(datetime([SA.Datenum],'convertfrom','datenum'));

% Get the quarters of the surveys
SurvQtr=1+floor(SurvMon/4);

% Get the winter(ONDJFM=1)/summer(AMJJAS=2) season of the surveys
SurvSsn=ones*size(SurvQtr);SurvSsn(SurvQtr == 2 | SurvQtr == 3)=2;

%fprintf('Stepping through Oct 1, (yr-1)-Sep 30,(yr) beach years...\n')

% step though beach years and get mean quarterly and annual profiles
byr=unique(BeachYr);
AnnZ=NaN(numel(byr),numel(SA(1).X1D));
AllQuarterlyZ=NaN(numel(byr)*4,numel(SA(1).X1D));
AllQuarterlyDatenum=NaN(1,numel(byr)*4);
AllQuarterlyBeachYears=NaN(1,numel(byr)*4);
AllQuarterlyQuarters=NaN(1,numel(byr)*4);
ny=0; % year counter
nq=0; % quarter counter

for y=byr
    ny=ny+1;
    
    % mean monthly profiles
    MonthlyZ=NaN(12,numel(SA(1).X1D));
    for m=1:12
        idx=find(BeachYr ==  y & SurvMon == m);
        if ~isempty(idx) 
            if numel([SA(idx).Z1Dtransect]) > 0
        zm=mean(vertcat(SA(idx).Z1Dtransect),1,'omitnan');
        MonthlyZ(m,:)=mean(vertcat(SA(idx).Z1Dtransect),1,'omitnan');
            end
        end
    end
    
    % mean quarterly profiles
    QuarterlyZ=NaN(4,numel(SA(1).X1D));
    for q=1:4
        n=1+3*(q-1);
        QuarterlyZ(q,:)=mean(MonthlyZ(n:n+2,:),1,'omitnan');
    end
    
    % Add QtrlyZ to All QtrlyZ sequential array of quarters over all years.
    %  Because these are beach years, the Oct-Dec 4th quarter
    %  is from the previous calendar year, so save first to
    %  make an all quarterly time series in the right order.
    for q=[4 1 2 3]
        nq=nq+1;
        AllQuarterlyZ(nq,:)=QuarterlyZ(q,:);
        AllQuarterlyBeachYears(nq)=y;
        AllQuarterlyQuarters(nq)=q;
        if q == 4
          AllQuarterlyDatenum(nq)=datenum(y-1,11,15); % middle OND qtr, previous year
        else
          AllQuarterlyDatenum(nq)=datenum(y,(q-1)*3+2,15); % middle qtrs
        end
    end
    
    % mean 2 season profiles
        SsnZ=NaN(2,numel(SA(1).X1D));
        SsnZ(1,:)=mean(QuarterlyZ([1 4],:),1,'omitnan'); % fall-winter
        SsnZ(2,:)=mean(QuarterlyZ([2 3],:),1,'omitnan'); % spring-summer
    
    % annual
        AnnZ(ny,:)=mean(SsnZ,1); % will = NaN if one or both seasons have no data
end
        
% global mean

GlobalZ=mean(AnnZ,1,'omitnan');

%idx=find(sum(~isnan(AnnZ)') > 0); % rows (beach years) with annual means
% BeachYrs=byr(idx);
% AnnZ=AnnZ(idx,:);

GBP(nmop).Mop=MopNumber;
GBP(nmop).X=SA(1).X1D;
GBP(nmop).QuarterlyDatenums=AllQuarterlyDatenum;
GBP(nmop).QuarterlyBeachYears=AllQuarterlyBeachYears;
GBP(nmop).QuarterlyQuarters=AllQuarterlyQuarters;
GBP(nmop).Zquarterly=AllQuarterlyZ;
GBP(nmop).AnnualBeachYears=byr;
GBP(nmop).Zannual=AnnZ;
GBP(nmop).Zglobal=GlobalZ;

end

%% remove bad global data seaward of msl at mop  513 
GBP(7).Zglobal(178:end)=NaN;

%% Syncronize the global profile x origins relative to
%   mop 510's MLW-MHW beach face location

zsmin=0.218;zsmax=1.344;% sync z range of mlw-mhw in navd88
% loop through mops and fit this global profile elev range to 
%  mop 510 by shifting x axis for each
for n=1:8
  [rmsd,xlag,xlength]=GetRmsdGMP(GBP,510,506+n,zsmin,zsmax);
  % rmsd = root mean square difference of fitted profiles for z range
  % xlag = shift in x in whole meters of mop to sync with 510
  % xlength = xshore length of overlap area
  fprintf('Mop %i  xlag= %i\n',GBP(n).Mop,xlag)
  % add global xshore sync info to struct array 
  GBP(n).Xlag=xlag;
  GBP(n).SyncMop=510;
  GBP(n).SyncZrange=[zsmin zsmax];
  GBP(n).Rmsd=rmsd;
end

%% save global mean profile struct array
matfile=['Mops' num2str(MopStart) 'to'  num2str(MopEnd) 'GBP.mat'];
fprintf('Saving GBP struct array in %s\n',matfile);
save(matfile,'GBP');

%% example plots using struct array fields
figure('position',[ 440    63   663   734]);
n=1; % choose a mop index 506+n = Mop number
MopName=['Mop: ' num2str(GBP(n).Mop)];
subplot(3,1,1);surf(GBP(n).X,GBP(n).QuarterlyDatenums,GBP(n).Zquarterly);
set(gca,'xdir','reverse');datetick('y');shading flat;BeachColorbar;
title([MopName ' Quarterly Means']);ylabel('Calendar Year');
view(2);xlabel('Xshore Distance (m)');

subplot(3,1,2);surf(GBP(n).X,GBP(n).AnnualBeachYears+.5,GBP(n).Zannual);
set(gca,'xdir','reverse');shading flat;BeachColorbar;
yl=[min(GBP(n).AnnualBeachYears) max(GBP(n).AnnualBeachYears)+1];
title([MopName ' Annual Means']);set(gca,'ylim',yl);
ylabel('Oct-Sep Beach Year');view(2);xlabel('Xshore Distance (m)')

subplot(3,1,3);plot(GBP(n).X,GBP(n).Zglobal,'k-','linewidth',2);
set(gca,'xdir','reverse');title([MopName ' Global Mean']);grid on;
xlabel('Xshore Distance (m)');ylabel('Elevation (m,NAVD88)');

function SA=EditSA(SA)

MopNumber=SA(1).Mopnum;

% remove bad transect profile jetski data from Mop 534 and 535
%   in the 7 feb 2005 jumbo survey
if MopNumber == 534 ||  MopNumber == 535
   idx=find([SA.Datenum] == datenum(2005,2,7));
   % remove bad deep section
   nbad=find( SA(idx).Z1Dtransect < -1.4);
   SA(idx).Z1Dtransect(nbad)=NaN;
end

% for mops 555 and 556, replace bad jet ski data where x < 250 and z < -1.7m on
%  4 Apr 2016 with linearly interpolated data

if MopNumber == 555 || MopNumber == 556
   idx=find([SA.Datenum] == datenum(2016,4,4));
   nbad=find(SA(idx).X1D < 250 & SA(idx).Z1Dtransect < -1.7);
   SA(idx).Z1Dtransect(nbad)=NaN;
   igood=find(~isnan(SA(idx).Z1Dtransect));
   SA(idx).Z1Dtransect=interp1(SA(idx).X1D(igood),...
       SA(idx).Z1Dtransect(igood),SA(idx).X1D);
end

end
